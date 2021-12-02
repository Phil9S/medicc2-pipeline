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

getMediccTableSeg <- function(object){
    valid = apply(object,MARGIN = 1,FUN = function(x) !any(is.na(x)))
    CN = object[valid,]
    CN = round(CN)
    #valid = binsToUseInternal(object)
    abs_valid <- abs_data[valid,]
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
        sample = colnames(object)[[sample_index]]
        cn_values = cn[start(rl)]
        chromosomeValues = chromosome[start(rl)]
        gr = GenomicRanges::GRanges(seqnames = chromosomeValues, IRanges(start(rl),
                                                                         width=runLength(rl), sample=sample, copynumber = cn_values))
        df = data.frame(gr) %>% dplyr::rename("chrom"="seqnames", "sample_id"="sample", "cn_a"=copynumber) %>% dplyr::select(-strand, -width) %>%
            dplyr::relocate(sample_id, chrom, start, end, cn_a)
        return(df)
    }
    dfs = dplyr::bind_rows(lapply(1:ncol(object), export_single_sample))
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

getBinsStartsEnds <- function(window=500000,chr,lengthChr){
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
seg_data <- read.table(args[1],header = TRUE,sep = "\t")

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

allchr=c(1:22) #Add 23 if want to include chrX
lengthChr <- c("1"="249250621","2"="243199373","3"="198022430",
               "4"="191154276","5"="180915260","6"="171115067",
               "7"="159138663","8"="146364022","9"="141213431",
               "10"="135534747","11"="135006516","12"="133851895",
               "13"="115169878","14"="107349540","15"="102531392",
               "16"="90354753","17"="81195210","18"="78077248",
               "19"="59128983","20"="63025520","21"="48129895",
               "22"="51304566")

posBins <- lapply(allchr,function(chr) 
    getBinsStartsEnds(window=3000000, chr, lengthChr[chr]))

if(type == "tot"){
    seg_data <- seg_data %>%
                    mutate(segVal = replace(segVal, segVal < 0, 0)) %>%
                    select(segcols)
    bins <- as.data.frame(getCNbins(posBins = posBins,data = seg_data,samples = unique(seg_data$sample)))
    bins <- data.frame((lapply(bins,FUN = function(x) as.numeric(x))),check.names = FALSE)
    
    processedbins <- do.call(rbind,lapply(posBins,FUN = function(x){
        chr <- x$chromosome
        start <- x$starts
        end <- x$ends
        data.frame(chromosome=chr,start=start,end=end,row.names = paste0(chr,":",start,"-",end))
    }))
    
    abs_data <- as.data.frame(cbind(processedbins,bins))
} else if(type == "as"){
    seg_dataA <- seg_dataAS %>%
                    select(all_of(segcolsAS)) %>%
                    select(-segValB) %>%
                    rename("segVal"="segValA")
    
    seg_dataB <- seg_dataAS %>%
                    select(all_of(segcolsAS)) %>%
                    select(-segValA) %>%
                    rename("segVal"="segValB")
    
    binsA <- as.data.frame(getCNbins(posBins = posBins,data = seg_dataA,samples = unique(seg_dataA$sample)))
    binsA <- data.frame((lapply(bins,FUN = function(x) as.numeric(x))),check.names = FALSE)
    
    binsB <- as.data.frame(getCNbins(posBins = posBins,data = seg_dataB,samples = unique(seg_dataB$sample)))
    binsB <- data.frame((lapply(bins,FUN = function(x) as.numeric(x))),check.names = FALSE)
    
    processedbins <- do.call(rbind,lapply(posBins,FUN = function(x){
        chr <- x$chromosome
        start <- x$starts
        end <- x$ends
        data.frame(chromosome=chr,start=start,end=end,row.names = paste0(chr,":",start,"-",end))
    }))
    
    abs_dataA <- as.data.frame(cbind(processedbins,binsA))
    abs_dataB <- as.data.frame(cbind(processedbins,binsB))
    
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
            object <- as.data.frame(abs_data[,colnames(abs_data) %in% y])
            table <- getMediccTableSeg(object = object)
            write.table(x = table,file = paste0(outfolder,patient_name),quote = F,sep = "\t",row.names = F,col.names = T)
        } else if(type == "as"){
            objectA <- as.data.frame(abs_dataA[,colnames(abs_dataA) %in% y])
            tableA <- getMediccTableSeg(object = objectA)
            objectB <- as.data.frame(abs_dataB[,colnames(abs_dataB) %in% y])
            tableB <- getMediccTableSeg(object = objectB)
            tableA$cn_b <- tableB$cn_a
            write.table(x = tableA,file = paste0(outfolder,patient_name),quote = F,sep = "\t",row.names = F,col.names = T)
        }
    }
})

## END
