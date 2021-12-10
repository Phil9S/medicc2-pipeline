# get a medicc compatible data frame from segment table
# Function getMediccTable() by Michael Schneider
# Modifed by Philip Smith
# Function getBinStartsEnds(),getCNbins(),getCNbins.bin(),getCNbins.sample() from CNpare package
# Modified by Philip Smith

args = commandArgs(trailingOnly=TRUE)
# Disable sci notation
options(scipen=999)

# Required packages
require(GenomicRanges, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(tidyr, quietly = TRUE, warn.conflicts = FALSE)

getMediccTableSeg <- function(object){
    valid = apply(object,MARGIN = 1,FUN = function(x) !any(is.na(x)))
    CN = object[valid,-which(colnames(object) %in% c("chromosome","start","end"))]
    CN = round(CN)
    #valid = binsToUseInternal(object)
    abs_valid <- object[valid,]
    info = rownames(abs_valid)
    #seqName = unlist(lapply(strsplit(info, ":"),[,1]))
    seqName = unlist(lapply(strsplit(info, ":"),FUN = function(x) x[1]))
    chrOrder = base::rle(seqName)[["values"]]
    chromosome <- factor(seqName, levels=chrOrder)
    ids = paste0("chr",chromosome,":", apply(CN, 1, function(x){
        v=paste(x, sep="-", collapse="-"); return(v)}))
    rl = Rle(ids[valid])
    export_single_sample <- function(sample_index){
        cn = CN[,sample_index]
        sample = colnames(CN)[[sample_index]]
        cn_values = cn[start(rl)]
        chromosomeValues = chromosome[start(rl)]
        gr = GenomicRanges::GRanges(seqnames = chromosomeValues, IRanges(start(rl),
                                                                         width=runLength(rl), sample=sample, copynumber = cn_values))
        df = data.frame(gr) %>% dplyr::rename("chrom"="seqnames", "sample_id"="sample", "cn_a"=copynumber) %>% dplyr::select(-strand, -width) %>%
            dplyr::relocate(sample_id, chrom, start, end, cn_a)
        return(df)
    }
    dfs = dplyr::bind_rows(lapply(1:ncol(CN), export_single_sample))
    # replace index with actual position
    dfs$start <- abs_valid$start[dfs$start]
    dfs$end <- abs_valid$end[dfs$end]
    return(dfs)
}

getCNbins<- function(posBins,data,samples){
    pb=data.table::rbindlist(posBins)
    out<-list()
    out<-lapply(seq_len(nrow(pb)), function(b) getCNbins.bin(b,pb,data,samples))
    cn<-as.matrix(do.call(rbind,out))
    colnames(cn)<-samples
    return(cn)
}

getBinsStartsEnds <- function(window,chr,lengthChr){
    divideChr <- seq(0, lengthChr, window)
    starts <- divideChr[-c(length(divideChr))] + 1
    ends <- divideChr[-c(1)]
    
    binsSe <- list(chromosome=chr,starts=starts,ends=ends)
    binsSe
}

getCNbins.bin <- function(b,pb,data,samples){
    out<-list()
    chrom<-as.character(pb[b,1])
    start<-as.numeric(pb[b,2])
    end<-as.numeric(pb[b,3])
    cn<-data[(data$chromosome%in%chrom & data$start<=start & data$end>=end),]
    out<-lapply(samples, function(s) getCNbins.sample(s,cn))
    out<-do.call(cbind,out)
    colnames(out)<-samples
    return(out)
}

getCNbins.sample <- function(s,cn){
    if (nrow(cn)!=0){
        segVal <- cn[cn$sample==s, "segVal"]
        segVal <- ifelse(length(segVal)!=0, segVal, NA)
    } else {
        segVal <- NA
    }
    return(segVal)
}

## Set id cols
metacols <- c("PATIENT_ID","SAMPLE_ID")
segcols <- c("chromosome","start","end","segVal","sample")
segcolsAS <- c("chromosome","start","end","segValA","segValB","sample")

# Load sample CN
seg_data <- read.table(args[1],header = TRUE,sep = "\t",stringsAsFactors = F)

if(any(grepl(pattern = "chr*",seg_data$chromosome))){
    seg_data$chromosome <- gsub(pattern = "chr",replacement = "",seg_data$chromosome)
}

if(any(grepl(pattern = "X",seg_data$chromosome))){
    seg_data$chromosome[seg_data$chromosome == "X"] <- 23
}

if(is.character(seg_data$chromosome)){
    chrs <- as.character(seq.int(1,23,1))
    seg_data <- seg_data[seg_data$chromosome %in% chrs,]
    seg_data$chromosome <- as.numeric(seg_data$chromosome)
}


# Set output folder
outfolder <- args[2]

# Load sample meta.data
meta.data <- read.table(args[3],header = TRUE,sep = "\t")
meta.data <- meta.data[,which(colnames(meta.data) %in% metacols)]

# Is allele-specific
if(ncol(seg_data) == 5){
    if(all(colnames(seg_data) %in% segcols)){
        type <- "tot"    
    } else {
        stop(paste0("Incorrectly named columns in segment table, column names should be ",paste0(segcols,collapse = ", ")))
    }
} else if(ncol(seg_data) == 6){
    if(all(colnames(seg_data) %in% segcolsAS)){
        type <- "as"    
    } else {
        stop(paste0("Incorrectly named columns in segment table, column names should be ",paste0(segcolsAS,collapse = ", ")))   
    }
} else {
    stop("Incorrect number of columns in segment table")
}

## set bin size
bin <- 30000

## Get chromosome lengths and names
lengthChr <- unlist(lapply(split(seg_data,f = seg_data$chromosome),FUN = function(x) max(x$end)))
lengthChr <- lengthChr[which(names(lengthChr) %in% unique(seg_data$chromosome))]
allchr <- names(lengthChr)

## Generate bin positions from chromsome lengths and bin size
posBins <- lapply(allchr,function(chr) 
    getBinsStartsEnds(window=bin, chr, lengthChr[chr]))

if(type == "tot"){
    seg_dataTOT <- seg_data %>%
                    mutate(segVal = replace(segVal, segVal < 0, 0)) %>%
                    select(all_of(segcols))
    
    bins <- as.data.frame(getCNbins(posBins = posBins,data = seg_dataTOT,samples = unique(seg_dataTOT$sample)))
    bins <- data.frame((lapply(bins,FUN = function(x) as.numeric(x))),check.names = FALSE)
    
    processedbins <- do.call(rbind,lapply(posBins,FUN = function(x){
        chr <- x$chromosome
        start <- x$starts
        end <- x$ends
        data.frame(chromosome=chr,start=start,end=end,row.names = paste0(chr,":",start,"-",end))
    }))
    
    abs_data <- as.data.frame(cbind(processedbins,bins))
} else if(type == "as"){
    # Generate new sample by allele table
    seg_dataA <- seg_data %>%
                    pivot_longer(cols = c(segValA,segValB)) %>%
                    mutate(sample = paste0(sample,"_",name)) %>%
                    arrange(sample,chromosome,start) %>%
                    select(-name) %>%
                    rename(segVal = "value") %>%
                    relocate(sample, .after = last_col())
    
    binsA <- as.data.frame(getCNbins(posBins = posBins,data = seg_dataA,samples = unique(seg_dataA$sample)))
    binsA <- data.frame((lapply(binsA,FUN = function(x) as.numeric(x))),check.names = FALSE)
    
    processedbins <- do.call(rbind,lapply(posBins,FUN = function(x){
        chr <- x$chromosome
        start <- x$starts
        end <- x$ends
        data.frame(chromosome=chr,start=start,end=end,row.names = paste0(chr,":",start,"-",end))
    }))
    
    abs_dataA <- as.data.frame(cbind(processedbins,binsA))
}


# Split sample lists by patient
sample_by_pat_list <- split(meta.data$SAMPLE_ID,f=as.factor(meta.data$PATIENT_ID))

# Perform iterative extraction of medicc compatible tables for each patient
lapply(names(sample_by_pat_list),FUN = function(x){
    # Extract patient 
    patient_name <- as.character(x)
    # vector of sample ids
    y <- unlist(sample_by_pat_list[[patient_name]])
    # if sample vector contains more than one sample then extract medicc
    # formatted table from segtable
    if(length(y) > 1){
        if(type == "tot"){
            object <- as.data.frame(abs_data[,colnames(abs_data) %in% c("chromosome","start","end",y)])
            table <- getMediccTableSeg(object = object)
            write.table(x = table,file = paste0(outfolder,patient_name),quote = F,sep = "\t",row.names = F,col.names = T)
        } else if(type == "as"){
            y2 <- paste0(y,"_",rep(c("segValA","segValB"),times=length(y)))
            objectA <- as.data.frame(abs_dataA[,colnames(abs_dataA) %in% c("chromosome","start","end",y2)])
            
            tableA <- getMediccTableSeg(object = objectA)
            
            tableAf <- tableA %>%
                        mutate(group = ifelse(grepl(sample_id,pattern = "_segValA"),"segValA","segValB")) %>%
                        mutate(sample_id = gsub(sample_id,pattern = "_segValA|_segValB",replacement = "")) %>%
                        group_by(sample_id) %>%
                        pivot_wider(names_from = group,values_from = cn_a) %>%
                        rename(cn_a = "segValA",cn_b = "segValB")

            write.table(x = tableAf,file = paste0(outfolder,patient_name),quote = F,sep = "\t",row.names = F,col.names = T)
        }
    }
})

## END
