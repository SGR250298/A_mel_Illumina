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
P450_results <- snakemake@output[["P450_results"]]

#P450 search
search_gff <- import.gff3(gff_file,feature.type = "CDS")
CDS_search <- as.data.table(mcols(search_gff))
searchP450 <- CDS_search[grepl("P450", product,ignore.case = TRUE) & !grepl("reductase", product),
           .(unique(protein_id))]
fwrite(searchP450, P450_results, col.names = FALSE)
sessionInfo()