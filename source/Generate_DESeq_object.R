#!/usr/bin/env Rscript

#catch log file in R
log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)
library(DESeq2)
library(rtracklayer)
library(tximport)

gff <- snakemake@input[["gff"]]
ddsfile <- snakemake@output[["dds"]]
quant_files <- snakemake@input[["quant_files"]]
sample_summary <- snakemake@input[["sample_summary"]]

#generate tx2gene approved format
mRNA <- import.gff3(gff, feature.type="mRNA")
gene_info <- as.data.table(mcols(mRNA))
tx2gene <- data.frame(gene_info[,.(Transcript=Name,Gene=as.character(gene))])

#import quant files and name
names(quant_files) <- basename(dirname(quant_files))
gene_counts <- tximport(quant_files,type="salmon", tx2gene = tx2gene)
#^excluded siRNAs and lncRNAs which were captured but not able to be mapped to transcripts (that's fine)

#read in sample data
sample_data <- fread(sample_summary)
sample_data[,hive_number:=paste0("hive",hive_number)]
#^put word in front of hive number to treat as characters
sample_data[,sample:=factor(sample,levels = names(quant_files))]
setorder(sample_data,sample)
#^match sample order to file listings
col_data <- data.frame(sample_data,row.names = "sample")

#generate DESeq object
dds <- DESeqDataSetFromTximport(gene_counts,col_data,~1)
saveRDS(dds,ddsfile)

#save session information (versions of packages etc.)
sessionInfo()