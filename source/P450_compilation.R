#!/usr/bin/env Rscript

#catch log file in R
log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")


library(rtracklayer)
library(data.table)

gff_file <- snakemake@input[["gff"]]
P450_proteinlist <- snakemake@input[["gene_list"]]
P450_compilation <- snakemake@output[["P450_list"]]
#gff_file <- "data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.gff"
#P450_proteinlist <- "output/P450_list/Amel_P450_proteins.txt"

#P450 search
search_gff <- import.gff3(gff_file,feature.type = "CDS")
CDS_search <- as.data.table(mcols(search_gff))
P450_proteins <- fread(P450_proteinlist, header = FALSE, col.names = "protein_id")
setkey(CDS_search, protein_id)
P450_genelist <- CDS_search[P450_proteins, unique(gene)]
writeLines(P450_genelist, P450_compilation)
sessionInfo()