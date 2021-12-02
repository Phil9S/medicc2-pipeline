# get a medicc compatible data frame from QDNAseq object (copynumber is based on copynumber slot)
# Function written by Michael Schneider
# Modified and debugged by Philip Smith

args = commandArgs(trailingOnly=TRUE)

# Required libraries
require(Biobase, quietly = TRUE, warn.conflicts = FALSE)
require(QDNAseqmod, quietly = TRUE, warn.conflicts = FALSE)
require(GenomicRanges, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)

getMediccTable <- function(object, value="copynumber"){
    CN = Biobase::assayDataElement(object, value)
    CN = round(CN)
    #valid = binsToUseInternal(object)
    valid = object@featureData@data$use
    info = rownames(object)
    #seqName = unlist(lapply(strsplit(info, ":"),[,1]))
    seqName = unlist(lapply(strsplit(info, ":"),FUN = function(x) x[1]))
    chrOrder = base::rle(seqName)[["values"]]
    chromosome <- factor(seqName, levels=chrOrder)
    ids = paste0("chr",chromosome,":", apply(CN, 1, function(x){
        v=paste(x, sep="-", collapse="-"); return(v)}))
    rl = Rle(ids[valid])
    export_single_sample <- function(sample_index){
        cn = CN[valid, sample_index]
        sample = colnames(object)[[sample_index]]
        cn_values = cn[start(rl)]
        chromosomeValues = chromosome[valid][start(rl)]
        gr = GenomicRanges::GRanges(seqnames = chromosomeValues, IRanges(start(rl),
                                                                         width=runLength(rl), sample=sample, copynumber = cn_values))
        df = data.frame(gr) %>% dplyr::rename("chrom"="seqnames", "sample_id"="sample", "cn_a"=copynumber) %>% dplyr::select(-strand, -width) %>%
            dplyr::relocate(sample_id, chrom, start, end, cn_a)
        # replace index with actual position
        df$start <- fData(object)$start[df$start]
        df$end <- fData(object)$start[df$end]
        return(df)
    }
    dfs = dplyr::bind_rows(lapply(1:ncol(object), export_single_sample))
    return(dfs)
}

# Load sample CN
#print(args)
abs_data <- readRDS(args[1])

# Set output folder
outfolder <- args[2]

# Load sample meta.data
meta.data <- read.table(args[3],header = TRUE,sep = "\t")

# Filter data and meta to match
abs_data <- abs_data[,which(colnames(abs_data) %in% meta.data$SAMPLE_ID)]

# Split sample lists by patient
sample_by_pat_list <- split(meta.data$SAMPLE_ID,f=as.factor(meta.data$PATIENT_ID))

# Perform iterative extraction of medicc compatible tables for each patient
medicc.table <- lapply(names(sample_by_pat_list),FUN = function(x){
    # Extract patient 
    patient_name <- as.character(x)
    # vector of sample ids
    y <- unlist(sample_by_pat_list[[patient_name]])
    # if sample vector contains more than one sample then extract medicc
    # formatted table from QDNAseq object
    if(length(y) > 1){
        object <- abs_data[,colnames(abs_data) %in% y]
        table <- getMediccTable(object = object,value = "segmented")
	      table$cn_a[table$cn_a < 0] <- 0
        write.table(x = table,file = paste0(outfolder,patient_name),quote = F,sep = "\t",row.names = F,col.names = T)
    }
    return(table)
})

## END
