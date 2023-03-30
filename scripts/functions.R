suppressMessages(library(GenomicRanges))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))

# Set options
options(dplyr.summarise.inform = FALSE)
options(warn=-1)

## Segment to bin format
get_binned_segment_table <- function(segtable = NULL,binsize = 30000,qdnaseq.format=FALSE){
  if(is.null(segtable)){
    stop("no segment table")
  }
  if(any(!colnames(segtable) %in% c("chromosome","start","end","segVal","sample"))){
      stop("bad column names")
  }
  # Get seq lengths
  seqlengths <- unlist(lapply(split(segtable,f=segtable$chromosome),function(x) max(x$end)))
  # Make Granges from data.frame format
  table_gr <- makeGRangesFromDataFrame(segtable,keep.extra.columns = T,ignore.strand = T)
  # Generate tiled genome of given bin size and seq.length
  tiles <- unlist(tileGenome(seqlengths = seqlengths,tilewidth = binsize))
  
  # Merge bins and segment tables to form annotated bins
  #binned <- as.data.frame(mergeByOverlaps(tiles,table_gr))
  binned <- mergeByOverlaps(tiles,table_gr)
  binned <- data.frame(binned$tiles,binned$table_gr,binned$segVal,binned$sample)
  
  cat(paste0(args[4],"[",format(Sys.time(),
    "%H:%M:%S"),"]","[MAIN] Binning samples across ",binsize/1000,"kb bins (this may take some time)\n"))

  binned_mat <- as.data.frame(binned %>%
                                  select(-c(contains(".1"),"strand","width","binned.segVal","binned.sample")) %>%
                                  distinct() %>%
                                  pivot_wider(values_from = "segVal",names_from = "sample",values_fn = median) %>%
                                  dplyr::rename(c("chromosome"="seqnames","start"="start","end"="end")))
  
  rowN <- paste0(binned_mat[,1],":",binned_mat[,2],"-",binned_mat[,3])
  binned_mat <- binned_mat[,-c(1:3)]
  binned_mat <- apply(binned_mat,MARGIN = 2,FUN = function(x) as.numeric(x))
  rownames(binned_mat) <- rowN
  naRows <- which(apply(binned_mat,1,FUN = function(x) any(is.na(x))))
  if(length(naRows) > 0){
    binned_mat <- binned_mat[-c(naRows),]  
  }
  if(!qdnaseq.format){
    cnames <- colnames(binned_mat)
    binned_mat <- as.data.frame(binned_mat)
    binned_mat <- cbind(str_split(rownames(binned_mat),pattern = ":|-",n = 3,simplify = T),binned_mat)
    rownames(binned_mat) <- NULL
    colnames(binned_mat) <- c("chromosome","start","end",cnames)
    binned_mat$chromosome <- factor(binned_mat$chromosome,levels = seq.int(1:23))
    binned_mat$start <- as.numeric(as.character(binned_mat$start))
    binned_mat$end <- as.numeric(as.character(binned_mat$end))
    return(binned_mat)
  } else {
    return(as.matrix(binned_mat))
  }
}

extract_qdnaseq <- function(data=NULL,value="segmented",qdnaseq.format=FALSE){
    if(is.null(data)){
        stop("no data")
    }
    if(!value %in% c("segmented","copynumber")){
        stop("unknown value")
    }
    if(class(data) != "QDNAseqCopyNumbers"){
        stop("QDNAseq class object not provided")
    }
    CN <- Biobase::assayDataElement(data, value)
    valid = data@featureData@data$use
    info = rownames(data)
    CN <- CN[valid,]
    if(!qdnaseq.format){
        cnames <- colnames(CN)
        CN <- as.data.frame(CN)
        CN <- cbind(str_split(rownames(CN),pattern = ":|-",n = 3,simplify = T),CN)
        rownames(CN) <- NULL
        colnames(CN) <- c("chromosome","start","end",cnames)
        CN$chromosome <- factor(CN$chromosome,levels = seq.int(1:23))
        CN$start <- as.numeric(as.character(CN$start))
        CN$end <- as.numeric(as.character(CN$end))
        return(CN)
    } else {
        return(as.matrix(CN))
    }
}

normalised_segments <- function(binned=NULL,mapping=NULL,type=NULL){
    if(is.null(binned)){
        stop("no binned data")
    }
    if(is.null(mapping)){
        stop("no mapping file")
    }
    sample_split <- split(mapping,f = mapping$PATIENT_ID)
    
    ## Set up parallel compute
    if(cores > 1){
	cat(paste0(args[4],"[",format(Sys.time(),"%H:%M:%S"),"]","[MAIN] Running parallel on ",cores," cores\n"))
        doParallel::registerDoParallel(cores=cores)
        collapsed_subsets <- foreach::foreach(x=sample_split) %dopar% {
            s <- subset_binned(binned,x,type)
        }
	stopImplicitCluster()
        names(collapsed_subsets) <- names(sample_split)
    } else {
        collapsed_subsets <- lapply(sample_split,FUN = function(x) subset_binned(binned,x,type))
    }
    return(collapsed_subsets)
}


subset_binned <- function(binned,x,type){
    cat(paste0(args[4],"[",format(Sys.time(),"%H:%M:%S"),"]","[MAIN] Normalising segments for ",unique(x$PATIENT_ID),"\n"))
    if(type == "segment_as"){
        samples <- as.character(outer(x$SAMPLE_ID,paste0(":",c("segValA","segValB")),paste0))
    } else {
        samples <- as.character(x$SAMPLE_ID)
    }
    #samples <- colnames(binned)[!colnames(binned) %in% c("chromosome","start","end")]
    binned_sub <- binned[,colnames(binned) %in% c("chromosome","start","end",samples)] %>%
    #binned_sub <- binned %>%    
        mutate(across(where(is.numeric), round)) %>%
        unite(col = "const",4:ncol(.))
    collapsed_segments <- collapse_bins(binned_sub,samples=samples)
    return(collapsed_segments)
}

collapse_bins <- function(binned_sub,samples=samples){
    
    binned_sub_collapse <- as.data.frame(binned_sub) %>%
        group_by(chromosome) %>%
        mutate(change = const != lag(const)) %>%
        mutate(change = ifelse(is.na(change),TRUE,change)) %>%
        ungroup() %>%
        mutate(segs = cumsum(change)) %>%
        group_by(chromosome,segs) %>%
        summarise(
            across(start,min),
            across(end,max),
            across(const,unique)) %>%
        ungroup() %>%
        separate(col = "const",into = samples,sep = "_") %>%
        dplyr::relocate(chromosome,start,end) %>%
        dplyr::select(-segs) %>%
        pivot_longer(cols = 4:ncol(.),names_to = "sample",values_to = "segVal") %>%
        arrange(sample) %>%
        dplyr::relocate(chromosome,start,end,segVal,sample) %>%
        mutate(segVal = as.numeric(segVal)) %>%
        distinct()
    
    # binned_sub_collapse <- as.data.table(binned_sub)[,as.data.table(reduce(IRanges(start, end))),
    #                                                  by = .(chromosome,const)] %>%
    #     separate(col = "const",into = samples,sep = "_") %>%
    #     select(-c(width)) %>%
    #     dplyr::relocate(chromosome,start,end) %>%
    #     pivot_longer(cols = 4:ncol(.),names_to = "sample",values_to = "segVal") %>%
    #     arrange(sample) %>%
    #     dplyr::relocate(chromosome,start,end,segVal,sample)
    return(binned_sub_collapse)
}

get_medicc_tables <- function(norm_segs=NULL,write=TRUE,outputdir=NULL,type=NULL){
    if(is.null(norm_segs)){
        stop("no segs provided")
    }
    if(is.null(outputdir)){
        stop("no output dir")
    }
    if(is.null(type)){
        stop("no type")
    }
    withdiploid <- lapply(norm_segs,FUN = function(x) medicc_format(x,type))
    if(write){
        lapply(names(withdiploid),FUN = function(x){
            file_name <- x
            filesave <- as.data.frame(withdiploid[[x]])
            #cat(paste0("write file ",file_name),"\n")
            write.table(filesave,paste0(outputdir,file_name),append = F,quote = F,sep = "\t",row.names = F,col.names = T)
        })
    }
    return(withdiploid)
}

medicc_format <- function(x,type){
    #y <- x[x$sample == x$sample[1],]
    #y$sample <- rep("diploid",nrow(y))
    #y$segVal <- as.numeric(rep(2,nrow(y)))
    #a <- rbind(x,y)
    a <- x 
    if(type == "segment_as"){
        colnames(a) <- c("chrom","start","end","cn_a","cn_b","sample_id")
        #a$chrom <- paste0("chr",a$chrom)
        a <- a[,c("sample_id","chrom","start","end","cn_a","cn_b")]
        a$cn_a[a$cn_a < 0] <- 0
        a$cn_b[a$cn_b < 0] <- 0
    } else {
        colnames(a) <- c("chrom","start","end","cn_a","sample_id")
        #a$chrom <- paste0("chr",a$chrom)
        a <- a[,c("sample_id","chrom","start","end","cn_a")]
        a$cn_a[a$cn_a < 0] <- 0
    }
    
    return(as.data.frame(a))
}

get_normalised_segments <- function(data=NULL,mapping=NULL,type=NULL,value="segmented",binsize=30000,qdnaseq.format=FALSE){
    if(is.null(data)){
        stop("no binned data")
    }
    if(is.null(mapping)){
        stop("no mapping file")
    }
    if(is.null(type)){
        stop("no type specfified")
    }
    if(!type %in% c("segment","segment_as","qdnaseq","snp")){
        stop("bad type")
    }
    switch(type,
    segment={
        binned_mat <- get_binned_segment_table(segtable = data,
                                               binsize=binsize,
                                               qdnaseq.format=qdnaseq.format)
        
        normalised_segs_tibble <- normalised_segments(binned = binned_mat,
                                                      mapping = mapping,
                                                      type=type)
    },
    segment_as={
        data.split <- split_alleles(data = data)
        binned_mat <- get_binned_segment_table(segtable = data.split,
                                               binsize=binsize,
                                               qdnaseq.format=qdnaseq.format)
        
        normalised_segs_tibble <- normalised_segments(binned = binned_mat,
                                                      mapping = mapping,
                                                      type=type)
        
        normalised_segs_tibble <- collapse_alleles(data = normalised_segs_tibble)
    },
    qdnaseq={
        binned_mat <- extract_qdnaseq(data = data,
                                      value = value,
                                      qdnaseq.format=qdnaseq.format)
        
        normalised_segs_tibble <- normalised_segments(binned = binned_mat,
                                                      mapping = mapping,
                                                      type=type)
    },
    snp={
        stop("do nothing")
    })
    
    normalised_segments <- as.list(normalised_segs_tibble)
    return(normalised_segments)
}

split_alleles <- function(data = NULL){
    if(is.null(data)){
        stop("no data")
    }
    data <- data %>%
        pivot_longer(cols = c("segValA","segValB"),values_to = "segVal") %>%
        mutate(sample = paste0(sample,":",name)) %>%
        relocate(sample,.after = last_col()) %>%
        select(-name) %>%
        arrange(sample,chromosome,start,end)
    return(data)
}

collapse_alleles <- function(data = NULL){
    if(is.null(data)){
        stop("no data")
    }
    collapse.data <- lapply(data,function(x){
        x <- x %>%
            separate(sample,sep = ":",into = c("sample","name")) %>%
            pivot_wider(names_from = name,values_from = segVal) %>%
            relocate(sample,.after = last_col()) %>%
            arrange(sample,chromosome,start,end)
        return(x)
            
    })
    return(collapse.data)
}

# Testing
## Load meta.data
# meta <- read.table("resources/mapping_file_example.tsv",header = T,sep = "\t")
# 
# ## Total seg input
# segs1 <- read.table("resources/segment_table_example.tsv",header = T,sep = "\t")
# test1 <- get_normalised_segments(data = segs1,mapping = meta,type = "segment")
# df1 <- test1[[2]]
# table(df1$sample)
# 
# ## QDNAseq input
# segs2 <- readRDS("resources/qdnaseqmod_example_file.rds")
# test2 <- get_normalised_segments(data = segs2,mapping = meta,type = "qdnaseq")
# df2 <- test2[[2]]
# table(df2$sample)
# 
# ## AS seg input
# segs3 <- read.table("resources/segment_table_allele_specfic_example.tsv",header = T,sep = "\t")
# test3 <- get_normalised_segments(data = segs3,mapping = meta,type = "segment_as")
# df3 <- test3[[2]]
# table(df3$sample)

# END
