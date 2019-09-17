library(DESeq2)
library(data.table)
library(rtracklayer)
library(ggplot2)
library(genefilter)
dds <- readRDS("output/040_DESeq/dds.Rds")
P450 <- readLines("output/P450_list/P450_compilation.txt")
sample_order <- c("H3_BF", "H2_OxFc", "H2_OxY", "H1b", "H2b", "H3b")
Set1 <- RColorBrewer::brewer.pal(9,"Set1")

#PCA plot
vst <- varianceStabilizingTransformation(dds,blind=TRUE)
plotPCA(vst,intgroup="exposure_status")

#filtering
row_mean <- rowMeans(counts(dds))
row_max <- rowMax(counts(dds))
dds_filtered <- dds[(row_mean>10 | row_max>20),
                    (dds$RIN>7 & dds$PCR_cycles<16 & dds$caste!="queen")]

#likelihood test
dds_lrt <- copy(dds_filtered)
dds_lrt$caste <- factor(dds_lrt$caste)
dds_lrt$exposure_status <- factor(dds_lrt$exposure_status)
design(dds_lrt) <- ~exposure_status + caste
dds_lrt <- DESeq(dds_lrt, test='LRT',reduced = ~1)
res_lrt <- results(dds_lrt,alpha=0.1)
save_lrt <- rownames(subset(res_lrt,!is.na(padj)))
vst_lrt <- varianceStabilizingTransformation(dds_lrt[save_lrt,],blind = FALSE)

#plotPCA(vst_lrt,intgroup="exposure_status")
vst_counts <- assay(vst_lrt)
vst_PCA <- prcomp(t(vst_counts),center = TRUE)
plot(vst_PCA$x)
percent_var <- ((vst_PCA$sdev)^2/(sum((vst_PCA$sdev)^2)))*100
vst_table <- data.table(percent_var,variable=paste0('PC',1:length(percent_var)))
PCA_dt <- data.table(vst_PCA$x,keep.rownames = TRUE)
PCA_plot <- melt(PCA_dt,id.vars = 'rn')
PCA_mod <- merge(PCA_plot,vst_table)
sample_table <- data.table(data.frame(colData(dds_lrt)),keep.rownames = TRUE)
plot_data <- merge(PCA_mod,sample_table,by='rn')
plot_data[,f_lab:=paste0(variable,' ',round(percent_var,1),'%')]
ggplot(plot_data[percent_var>5],aes(x=exposure_status,y=value,colour=hive_location))+geom_point()+facet_wrap(~f_lab)

#DE test
dds_treatment <- copy(dds_filtered)
dds_treatment$exposure_status <- factor(dds_treatment$exposure_status)
dds_treatment$caste <- factor(dds_treatment$caste)
design(dds_treatment) <- ~exposure_status
dds_treatment <- DESeq(dds_treatment)
plotMA(dds_treatment)
res <- results(dds_treatment,contrast = c("exposure_status","pyrethroid","organic"),lfcThreshold = 1,alpha = 0.1)
subset(res,padj<0.1 & rownames(res) %in% P450)
plotCounts(dds_treatment,"LOC551197",intgroup ="exposure_status")

#P450 heatmap
topVarGenes <- head(order(rowVars(assay(vst_lrt)), decreasing = TRUE), 50)
sig_genes <- rownames(res_lrt[order(res_lrt$padj),][1:50,])
mat  <- assay(vst_lrt)[ rownames(vst_lrt) %in% P450, ]
mat  <- t(scale(t(mat),center = TRUE))
hc <- hclust(dist(mat,method = "minkowski"),method = "ward.D2")
gene_order <- hc$labels[hc$order]
anno <- as.data.frame(colData(vst_lrt)[, c("caste","exposure_status")])
anno_dt <- data.table(anno, keep.rownames = TRUE)
anno_dt[,rn:= factor(rn, levels = sample_order)]
anno_dt[,Xcol:= Set1[as.numeric(as.factor(exposure_status))]]
setkey(anno_dt, rn)
mat_long <- as.data.table(melt(mat))
mat_long[,Var2:= factor(as.character(Var2),levels=sample_order)]
mat_long[,Var1:= factor(as.character(Var1),levels=gene_order)]
P450_heatmap <- ggplot(mat_long, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster() + 
  scale_fill_viridis_c(guide = guide_colourbar(title = "Transformed\nRead Counts")) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme_minimal(base_size = 10, base_family = "Times") +
  theme(axis.text.x = element_text(colour = anno_dt[,Xcol]),axis.text.y = element_text(face = "italic"))
ggsave("P450.pdf",P450_heatmap,width = 150, height = 150, units = "mm")

#General Heatmap
mat2  <- assay(vst_lrt)[ rownames(vst_lrt) ]
mat2  <- t(scale(t(mat2),center = TRUE))
hc2 <- hclust(dist(mat2,method = "minkowski"),method = "ward.D2")
gene_order2 <- hc2$labels[hc2$order]
anno2 <- as.data.frame(colData(vst_lrt)[, c("caste","exposure_status")])
anno_dt2 <- data.table(anno2, keep.rownames = TRUE)
anno_dt2[,rn:= factor(rn, levels = sample_order)]
anno_dt2[,Xcol:= Set1[as.numeric(as.factor(exposure_status))]]
setkey(anno_dt2, rn)
mat_long2 <- as.data.table(melt(mat2))
mat_long2[,Var2:= factor(as.character(Var2),levels=sample_order)]
mat_long2[,Var1:= factor(as.character(Var1),levels=gene_order)]
General_heatmap <- ggplot(mat_long2, aes(x=Var2, y=Var1, fill=value)) + 
  geom_raster() + 
  scale_fill_viridis_c(guide = guide_colourbar(title = "Transformed\nRead Counts")) + 
  xlab(NULL) + 
  ylab(NULL) +
  theme_minimal(base_size = 10, base_family = "Times") +
  theme(axis.text.x = element_text(colour = anno_dt[,Xcol]),axis.text.y = element_text(face = "italic"))
ggsave("General.pdf",P450_heatmap,width = 150, height = 150, units = "mm")