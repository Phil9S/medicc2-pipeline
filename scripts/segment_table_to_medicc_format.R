# get a medicc compatible data frame from segment table
args = commandArgs(trailingOnly=TRUE)
# Disable sci notation
options(scipen=999)

# Required packages
require(GenomicRanges, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(tidyr, quietly = TRUE, warn.conflicts = FALSE)
require(stringr, quietly = TRUE, warn.conflicts = FALSE)

source("scripts/functions.R")

# Set id cols
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

# set bin size
bin <- 30000

# Set output folder
outfolder <- args[2]

# Load sample meta.data
meta.data <- read.table(args[3],header = TRUE,sep = "\t")
meta.data <- meta.data[,which(colnames(meta.data) %in% metacols)]

norm_segs <- get_normalised_segments(data = seg_data,
		binsize=bin,
		mapping = meta.data,
		type = "segment")

medicc_out <- get_medicc_tables(norm_segs = norm_segs,outputdir=outfolder,write=TRUE)

outdir_up <- gsub("input_files/","",outfolder)
saveRDS(norm_segs,paste0(outdir_up,"norm_segs.rds"))
saveRDS(medicc_out,paste0(outdir_up,"medicc_out.rds"))

## END
