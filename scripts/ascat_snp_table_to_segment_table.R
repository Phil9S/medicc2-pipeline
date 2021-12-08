# get a allele normalised segment table from ASCAT SNP input table
# Written by Philip Smith

args = commandArgs(trailingOnly=TRUE)

# Disable sci notation
options(scipen=999)

require(GenomicRanges, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)

#getSegTableAscatRs
getSegTableAscatRs <- function(object){
    CN = object[,c(4:5)]
    CN = round(CN)
    info = object[,c(1:3)]
    seqName = object[,1]
    chrOrder = base::rle(seqName)[["values"]]
    chromosome <- factor(seqName, levels=chrOrder)
    ids = paste0("chr",chromosome,":", paste0(CN[,1],"-",CN[,2]))
    rl = Rle(ids)
    export_single_sample <- function(sample_index){
        cn = CN[,sample_index]
        sample = colnames(CN)[[sample_index]]
        cn_values = cn[start(rl)]
        chromosomeValues = chromosome[start(rl)]
        gr = GenomicRanges::GRanges(seqnames = chromosomeValues,
                                    IRanges(start(rl),
                                    width=runLength(rl),
                                    sample=sample,
                                    copynumber = cn_values))
        df = data.frame(gr) %>%
            dplyr::rename("chrom"="seqnames", "sample_id"="sample", "cn_a"=copynumber) %>%
            dplyr::select(-strand, -width) %>%
            dplyr::relocate(sample_id, chrom, start, end, cn_a)
        # replace index with actual position
        df$start <- object$start[df$start]
        df$end <- object$start[df$end]
        return(df)
    }
    dfs = dplyr::bind_rows(lapply(1:2, export_single_sample))
    dfsA <- dfs[dfs$sample_id == "cn_a",]
    dfsB <- dfs[dfs$sample_id == "cn_b",]
    dfsA$cn_b <- dfsB$cn_a
    dfs <- dfsA %>%
            select(-sample_id)
    return(dfs)
}

filename <- gsub(pattern = "\\..*",replacement="",basename(args[1]))
t <- read.table(args[1],header = T,sep = "\t")
c <- t %>%
        select(Chromosome,Position,Copy.number,Minor.allele) %>%
        mutate(cn_a = Copy.number - Minor.allele) %>%
        mutate(end = Position) %>%
        select(-Copy.number) %>%
        select(Chromosome,Position,end,cn_a,Minor.allele) %>%
        rename(start="Position",cn_b="Minor.allele")
        
table <- getSegTableAscatRs(c)
write.table(table,file = paste0(filename,".segments.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)

## END
