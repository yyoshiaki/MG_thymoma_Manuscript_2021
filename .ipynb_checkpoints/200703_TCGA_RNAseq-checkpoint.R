# 2020/07/03 Yoshiaki Yasumizu
# https://rpubs.com/tiagochst/TCGAbiolinks_to_DESEq2
# https://support.bioconductor.org/p/87333/

library(TCGAbiolinks)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(EnhancedVolcano)
library(viridis)
library(plyr)
library(vsn)
library(pals)
library(tidyverse)
library(dplyr)
library(scales)
library(grid)

setwd("/mnt/media32TB/home/yyasumizu/bioinformatics/TCGA_thymoma")
theme_set(theme_classic(base_size = 18, base_family = "Helvetica"))

proj <- "TCGA-THYM"
query <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# row.names(data)[1:10]
# row.names(data) <- getBM(filters="ensembl_gene_id", attributes="hgnc_symbol", values=row.names(data), mart=mart)$hgnc_symbol
# getBM(filters="ensembl_gene_id", attributes="hgnc_symbol", values=row.names(data)[1:10], mart=mart)$hgnc_symbol
row.names(data) <- rowRanges(data)$external_gene_name

# both are normal portion. remove theses.
colnames(data[,duplicated(substr(colnames(data), 1, 12))])

data <- data[,!duplicated(substr(colnames(data), 1, 12))]
colnames(data) <- substr(colnames(data), 1, 12)
# data[,duplicated(colnames(data))]

write.csv(assay(data), '201123_TCGATHYM_HTSeq_raocounts.csv')

# this is only indexed data.
clinical <- GDCquery_clinic(project = "TCGA-THYM", type = "clinical")
query_cli <- GDCquery(project = "TCGA-THYM", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query_cli)
clinical <- GDCprepare_clinic(query_cli, clinical.info = "patient")
clinical <- clinical[!duplicated(clinical$bcr_patient_barcode),]
rownames(clinical) <- clinical$bcr_patient_barcode
clinical <- clinical[data$bcr_patient_barcode,]

colData(data) <- cbind(colData(data), clinical[row.names(colData(data)),])
data <- data[!is.na(row.names(data)),!is.na(data$primary_pathology_history_myasthenia_gravis)]
data <- data[!duplicated(row.names(data)),]
data$primary_pathology_history_myasthenia_gravis <- 
  factor(data$primary_pathology_history_myasthenia_gravis, levels = c("NO", "YES"))

# examine the concordance of MG patients.
# clinical[clinical$primary_pathology_history_myasthenia_gravis=="YES",'bcr_patient_barcode']
# colnames(data[,data$primary_pathology_history_myasthenia_gravis == "YES"])

# ddsSE <- DESeqDataSet(data, design = ~ primary_pathology_history_myasthenia_gravis + 
#                         primary_pathology_histological_type_list)
ddsSE <- DESeqDataSet(data, design = ~ primary_pathology_history_myasthenia_gravis)
keep <- rowMeans(counts(ddsSE)) >= 5
ddsSE <- ddsSE[keep,]
ddsSE <- DESeq(ddsSE)

resultsNames(ddsSE)

res <- results(ddsSE, name = "primary_pathology_history_myasthenia_gravis_YES_vs_NO", alpha = 0.1, lfcThreshold = 1)
resOrdered <- res[order(res$padj),]

dea <- as.data.frame(resOrdered)
write.csv(dea, file = "200703_DSeq2WGCNA/DESeq2.res.csv")
summary(res)
dea
dea[c("NEFM", "NEFL", "NEFH", 'RYR3', 'GABRA5', 'GABRE', 'CHRNA1', 'RBFOX3'),]
plotMA(res, alpha = 0.1)

degs <- subset(resOrdered, padj < 0.1)
degs <- degs[abs(degs$log2FoldChange) > 1,]
genes.up <- row.names(degs[degs$log2FoldChange >= 0,])
genes.down <- row.names(degs[degs$log2FoldChange < 0,])
write(genes.up, file = "200703_DSeq2WGCNA/DESeq2.padj0.1.lfc1.upgenes.txt")
write(genes.down, file = "200703_DSeq2WGCNA/DESeq2.padj0.1lfc1.downgenes.txt")

# degs <- subset(resOrdered, padj < 0.05)
# degs <- degs[abs(degs$log2FoldChange) > 1,]
# genes.up <- row.names(degs[degs$log2FoldChange >= 0,])
# genes.down <- row.names(degs[degs$log2FoldChange < 0,])
# write(genes.up, file = "200703_DSeq2WGCNA/DESeq2.padj0.05lfc1.upgenes.txt")
# write(genes.down, file = "200703_DSeq2WGCNA/DESeq2.padj0.05lfc1.downgenes.txt")

# pheatmap
vsd <- vst(ddsSE)
meanSdPlot(assay(vsd))
df <- as.data.frame(colData(ddsSE)[,c("primary_pathology_history_myasthenia_gravis", "primary_pathology_histological_type_list")])
row.names(df) <- colnames(ddsSE)
levels(df$primary_pathology_histological_type_list) <- c( "Thymoma; Type A", "Thymoma; Type AThymoma; Type AB", "Thymoma; Type AB", "Thymoma; Type B1", "Thymoma; Type B1Thymoma; Type B2",
                                                          "Thymoma; Type B2", "Thymoma; Type B2Thymoma; Type B3", "Thymoma; Type B3", "Thymoma; Type C" )
df$primary_pathology_histological_type_list <- revalue(df$primary_pathology_histological_type_list, c("Thymoma; Type A" = "A", "Thymoma; Type AThymoma; Type AB" = "A;AB", "Thymoma; Type AB" = "AB", "Thymoma; Type B1" = "B1", "Thymoma; Type B1Thymoma; Type B2" = "B1;B2",
                                                       "Thymoma; Type B2" = "B2", "Thymoma; Type B2Thymoma; Type B3" = "B2;B3", "Thymoma; Type B3" = "B3", "Thymoma; Type C" = "C"))

# Specify colors
# cpal <- viridis::plasma(9)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cpal <- gg_color_hue(9)
ann_colors = list(
  primary_pathology_history_myasthenia_gravis = c(YES = "red", NO = "white"),
  primary_pathology_histological_type_list = c(A = cpal[1], "A;AB"=cpal[2], AB=cpal[3], B1=cpal[4], "B1;B2"=cpal[5], B2=cpal[6], "B2;B3"=cpal[7], B3=cpal[8], C=cpal[9])
)
zvsd <- (assay(vsd) - rowMeans(assay(vsd))) / rowSds(assay(vsd))

pdf("200703_DSeq2WGCNA/plots/pheatmap.DESeq2.padj0.1.lfc1.pdf", width = 16, height = 28)
pheatmap(assay(vsd)[row.names(degs),], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")
dev.off()

pdf("200703_DSeq2WGCNA/plots/pheatmap.vsd.CHRN.pdf", width = 16, height = 8)
pheatmap(assay(vsd)[read.csv("./data/HGNC_group-173_CHRN.csv", stringsAsFactors = FALSE)$Approved.symbol,], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")
dev.off()

# pheatmap(zvsd[c("NEFM", "NEFL", "GABRA5", "GABRE", "RYR3", "CHRNA1", "RBFOX1", "RBFOX3", "GLRA4", "NGB",
#                       "PLXNB3"),], annotation_col=df, color = brewer.brbg(10), annotation_colors = ann_colors,
#          clustering_distance_rows = "correlation", clustering_method = "ward.D2")
pheatmap(assay(vsd)[c("NEFM", "NEFL", "GABRA5", "GABRE", "RYR3", "CHRNA1", "RBFOX1", "RBFOX3", "GLRA4", "NGB",
                      "PLXNB3", "LHFPL4"),], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")

# PCA
levels(vsd$primary_pathology_histological_type_list) <- c( "Thymoma; Type A", "Thymoma; Type AThymoma; Type AB", "Thymoma; Type AB", "Thymoma; Type B1", "Thymoma; Type B1Thymoma; Type B2",
                                                          "Thymoma; Type B2", "Thymoma; Type B2Thymoma; Type B3", "Thymoma; Type B3", "Thymoma; Type C" )
vsd$primary_pathology_histological_type_list <- revalue(vsd$primary_pathology_histological_type_list, c("Thymoma; Type A" = "A", "Thymoma; Type AThymoma; Type AB" = "A;AB", "Thymoma; Type AB" = "AB", "Thymoma; Type B1" = "B1", "Thymoma; Type B1Thymoma; Type B2" = "B1;B2",
                                                                                                      "Thymoma; Type B2" = "B2", "Thymoma; Type B2Thymoma; Type B3" = "B2;B3", "Thymoma; Type B3" = "B3", "Thymoma; Type C" = "C"))

pdf("200703_DSeq2WGCNA/plots/PCA.pdf", width = 10, height = 5)
p1 <- plotPCA(vsd, intgroup=c("primary_pathology_history_myasthenia_gravis")) + 
  stat_ellipse(type = "norm", level = 0.67, geom = "polygon", alpha = 0.1, aes(fill = primary_pathology_history_myasthenia_gravis)) +
  ggtitle("Myasthenia gravis") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
p2 <- plotPCA(vsd, intgroup=c("primary_pathology_histological_type_list")) +  
  ggtitle("Histological type") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
gridExtra::grid.arrange(p1, p2, nrow = 1)
dev.off()

pca <- plotPCA(vsd, intgroup=c("primary_pathology_history_myasthenia_gravis"), returnData = TRUE) 
pca$primary_pathology_history_myasthenia_gravis <- revalue(pca$primary_pathology_history_myasthenia_gravis, c("NO" = "nonMG", "YES" = "MG"))
pca.mean <- dplyr::as_tibble(pca) %>% 
  group_by(primary_pathology_history_myasthenia_gravis) %>%
  summarise(PC1mean = mean(PC1), PC2mean = mean(PC2))
# plotPCA(vsd, intgroup=c("primary_pathology_history_myasthenia_gravis")) +
gg.mg <- ggplot(pca, aes(x = PC1, y = PC2, color = primary_pathology_history_myasthenia_gravis, label = primary_pathology_history_myasthenia_gravis)) +
  geom_point() +
  stat_ellipse(type = "norm", level = 0.67, geom = "polygon", alpha = 0.1, aes(fill = primary_pathology_history_myasthenia_gravis)) +
  ggtitle("Myasthenia gravis") +
  geom_label(data = pca.mean, aes(x = PC1mean, y = PC2mean), show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.title=element_blank()) + 
  guides(col = guide_legend(nrow =  1))
gg.mg
ggsave(filename = "200703_DSeq2WGCNA/plots/PCA.MG.pdf", width = 5, height = 5)

pca <- plotPCA(vsd, intgroup=c("primary_pathology_histological_type_list"), returnData = TRUE) 
pca.mean <- dplyr::as_tibble(pca) %>% 
  group_by(primary_pathology_histological_type_list) %>%
  summarise(PC1mean = mean(PC1), PC2mean = mean(PC2))
gg.who <- ggplot(pca, aes(x = PC1, y = PC2, color = primary_pathology_histological_type_list, label = primary_pathology_histological_type_list)) +
  geom_point() +
  stat_ellipse(type = "norm", level = 0.67, geom = "polygon", alpha = 0.1, aes(fill = primary_pathology_histological_type_list)) +
  ggtitle("Histological type") +
  geom_label(data = pca.mean, aes(x = PC1mean, y = PC2mean), show.legend = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.title=element_blank()) + 
  guides(col = guide_legend(nrow =  1))
gg.who
ggsave(filename = "200703_DSeq2WGCNA/plots/PCA.WHO.pdf", width = 5, height = 5)

# gridExtra::grid.arrange(gg.mg, gg.who, ncol = 2)

# volcano plot
res.vol <- res[res$pvalue < 0.05,]

# all DEGs
pdf("200703_DSeq2WGCNA/plots/vol.alls.pdf", width = 20, height = 12)
res1 <- lfcShrink(ddsSE, contrast = c("primary_pathology_history_myasthenia_gravis", 'YES', 'NO'), 
                  res=res, type = 'normal')
EnhancedVolcano(res.vol,
                lab = rownames(res.vol),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 1,
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                cutoffLineWidth = 0,
                ylim = c(0,13),
                xlim = c(-20,10),
                col=c("black", "black", "black", "red3"),)
dev.off()

pdf("200703_DSeq2WGCNA/plots/vol.alls.annot30.pdf", width = 6, height = 6)
EnhancedVolcano(res.vol,
                lab = rownames(res.vol),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 1,
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                cutoffLineWidth = 0,
                selectLab = row.names(resOrdered[1:30,]),
                ylim = c(0,13),
                xlim = c(-10,10),
                col=c("black", "black", "black", "red3"),)
dev.off()

# selected genes
plot.genes <- c("NEFM", 'GABRA5', 'RYR3', "GABRE", "TRAJ24", "CHRNA1", "IL13R2")
pdf("200703_DSeq2WGCNA/plots/vol.selectedgenes.pdf", width = 6, height = 8)
EnhancedVolcano(res.vol,
                lab = rownames(res.vol),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 0,
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                selectLab = plot.genes,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                cutoffLineWidth = 0,
                col=c("black", "black", "black", "red3"),)
dev.off()

# stacked bar plot
mg.who <- df
mg.who <- tibble::as_tibble(mg.who)
# mg.who$primary_pathology_history_myasthenia_gravis
# mg.who <- mg.who %>% mutate(count = 1)
mg.who <- mg.who %>% 
  dplyr::group_by(primary_pathology_histological_type_list, primary_pathology_history_myasthenia_gravis) %>%
  dplyr::tally()　# %>%
  # dplyr::filter(n > 5)
  
ggplot(mg.who, aes(x = primary_pathology_histological_type_list, y = n, fill = primary_pathology_history_myasthenia_gravis)) + 
  geom_bar( stat='identity')
ggsave("200703_DSeq2WGCNA/plots/MG.WHO.bar.pdf", width = 7, height = 2)

# violin
g <- "NEFM"
plot.violin <- function(g) {
  d <- plotCounts(ddsSE, gene=g, intgroup="primary_pathology_history_myasthenia_gravis", 
                  returnData=TRUE)
  
  levels(d$primary_pathology_history_myasthenia_gravis) <- c( "NO", "YES" )
  d$primary_pathology_history_myasthenia_gravis <- revalue(d$primary_pathology_history_myasthenia_gravis, c("NO" = "MG (-)", "YES" = "MG (+)"))
  
  pval <- res[g,"pvalue"]
  
  p <- ggplot(d, aes(x=primary_pathology_history_myasthenia_gravis, y=count)) + 
    geom_violin(adjust = 1, aes(fill = primary_pathology_history_myasthenia_gravis)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs("") + 
    ggtitle(g) +
    xlab("") + 
    ylab("") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none', plot.margin = margin(.1,.1,.1,.1, "cm")) +
    ggplot2::annotate("text", x = 1.2, y = 1, label = paste("p = ", format.pval(pval)))+
    geom_boxplot(width=0.1, fill="white")
  
  #     ylab(bquote(~Log[10]~count)) +
  return(p)
}
  
# p <- plot.violin(g)
# p + ggplot2::annotate("text", x = 2.2, y = 1, label = paste("p = ", format.pval(pval)))

p <- gridExtra::grid.arrange(plot.violin("NEFM"), plot.violin("RYR3"), plot.violin("GABRA5"),
                        plot.violin("GRIN2A"), plot.violin("CHRNA1"), plot.violin("KCNC1"),
                        plot.violin("PLXNB3"), plot.violin("RBFOX3"), plot.violin("GLRA4"), nrow = 3)

ggsave("200703_DSeq2WGCNA/plots/violin.MG.genes.pdf", width = 7, height = 7, p)

# WGCNA
library(WGCNA)

rowSds(assay(vsd))
plot(rowMeans((assay(vsd)), rowSds(assay(vsd))))

n.top <- 3000
keep <- rank(-rowSds(assay(vsd))) < n.top
# "RBFOX3" %in% row.names(vsd[keep,])
oed <- vsd[keep,]

gene.names=rownames(oed)
trans.oed=t(assay(oed))
dim(trans.oed)
datExpr=trans.oed

datTraits <- data.frame(vsd$primary_pathology_history_myasthenia_gravis, vsd$primary_pathology_histological_type_list,
           vsd$gender, (vsd$age_at_diagnosis / 365))
row.names(datTraits) <- colnames(vsd)
colnames(datTraits) <- c("history_myasthenia_gravis", "histological_type",  
                         "gender", "age_at_diagnosis" )

datTraits <- binarizeCategoricalColumns(datTraits, dropFirstLevelVsAll = FALSE)
colnames(datTraits)
datTraits <- datTraits[,c("history_myasthenia_gravis.YES.vs.all", "histological_type.A.vs.all", "histological_type.A;AB.vs.all", "histological_type.B1.vs.all", "histological_type.B2.vs.all", "histological_type.B2;B3.vs.all",
            "histological_type.B3.vs.all", "histological_type.C.vs.all", "gender.male.vs.all", "age_at_diagnosis" )]

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
pdf("200703_DSeq2WGCNA/plots/WGCNA/traits.samples.pdf", width = 16, height = 8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf("200703_DSeq2WGCNA/plots/WGCNA/topology.pdf", width = 10, height = 5)
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=cex1, col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = 5;
adjacency = adjacency(datExpr, power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
table(mergedColors)

# sizeGrWindow(12, 9)
pdf("200703_DSeq2WGCNA/plots/WGCNA/dendromerge.pdf", width = 10, height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "200703_DSeq2WGCNA/02-networkConstruction-stepByStep.RData")


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# sizeGrWindow(10,6)
pdf("200703_DSeq2WGCNA/plots/WGCNA/heat.trait.pdf", width = 8, height = 5)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

nSelect = 800
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;

pdf("200703_DSeq2WGCNA/plots/WGCNA/TOMplot.png", width = 6, height = 6)
TOMplot(1-plotDiss, selectTree, selectColors, col=viridis(10))
dev.off()

# Define variable mg containing the mg column of datTrait
mg = as.data.frame(datTraits$history_myasthenia_gravis.YES.vs.all)
names(mg) = "myasthenia_gravis"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, mg, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(mg), sep="")
names(GSPvalue) = paste("p.GS.", names(mg), sep="")

module.mg = "yellow"
column = match(module.mg, modNames);
moduleGenes = moduleColors==module.mg;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module.mg, "module"),
                   ylab = "Gene significance for myasthenia gravis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module.mg)

write(gene.names[moduleColors==module.mg], file = paste("200703_DSeq2WGCNA/",module.mg,".txt"), sep = "")

hist(res[gene.names[moduleColors==module.mg],'log2FoldChange'], breaks = 20)

write(gene.names[moduleColors==module.mg][res[gene.names[moduleColors==module.mg],'log2FoldChange'] > 1],
      file = paste("200703_DSeq2WGCNA/",module.mg,"lfc1.txt"), sep = "")

names(table(mergedColors))

for (c in names(table(mergedColors))) {
  column = match(c, modNames);
  moduleGenes = moduleColors==c
  write(gene.names[moduleColors==c], file = paste("200703_DSeq2WGCNA/",c,".txt", sep = ""))
}


library(magrittr)
library(clusterProfiler)

# reactome
library(ReactomePA)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <- gene.names[moduleColors==module.mg][res[gene.names[moduleColors==module.mg],'log2FoldChange'] > 1]
r <- res[genes,]
q <- getBM(filters="hgnc_symbol", attributes=c("entrezgene_id", "hgnc_symbol"), values=row.names(r), mart=mart, )
genes.entrez <- q$entrezgene_id
r <- r[q$hgnc_symbol,]
r <- r[!is.na(q$entrezgene_id),]
row.names(r) <- q[!is.na(q$entrezgene_id),'entrezgene_id']
rlfc <- r$log2FoldChange
names(rlfc) <- q[!is.na(q$entrezgene_id),'entrezgene_id']

x <- enrichPathway(gene=genes.entrez ,pvalueCutoff=0.2, readable=T)
head(as.data.frame(x))

barplot(x, showCategory=6)

pdf("200703_DSeq2WGCNA/plots/WGCNA/yellow.lfc1.reactome.dot.pdf", width = 10, height = 4)
dotplot(x, showCategory=6)
dev.off()

emapplot(x)

pdf("200703_DSeq2WGCNA/plots/WGCNA/yellow.lfc1.reactome.cnet.pdf", width = 8, height = 5)
cnetplot(x, showCategory = 6, categorySize="pvalue", foldChange=rlfc,)
dev.off()

genes <- gene.names[moduleColors==module.mg][res[gene.names[moduleColors==module.mg],'log2FoldChange'] < -1]
r <- res[genes,]
q <- getBM(filters="hgnc_symbol", attributes=c("entrezgene_id", "hgnc_symbol"), values=row.names(r), mart=mart, )
genes.entrez <- q$entrezgene_id
r <- r[q$hgnc_symbol,]
r <- r[!is.na(q$entrezgene_id),]
row.names(r) <- q[!is.na(q$entrezgene_id),'entrezgene_id']
rlfc <- r$log2FoldChange
names(rlfc) <- q[!is.na(q$entrezgene_id),'entrezgene_id']

x <- enrichPathway(gene=genes.entrez ,pvalueCutoff=0.2, readable=T)
head(as.data.frame(x))

barplot(x, showCategory=6)

pdf("200703_DSeq2WGCNA/plots/WGCNA/yellow.lfc1.down.reactome.dot.pdf", width = 10, height = 4)
dotplot(x, showCategory=6)
dev.off()

emapplot(x)

pdf("200703_DSeq2WGCNA/plots/WGCNA/yellow.lfc1.down.reactome.cnet.pdf", width = 8, height = 5)
cnetplot(x, categorySize="pvalue", foldChange=rlfc,)
dev.off()


for (c in names(table(mergedColors))) {
  genes <- read.csv(file = paste("200703_DSeq2WGCNA/",c,".txt", sep = "", header=FALSE))
  q <- getBM(filters="hgnc_symbol", attributes=c("entrezgene_id", "hgnc_symbol"), values=genes, mart=mart, )
  genes.entrez <- q$entrezgene_id
  genes.entrez <- genes.entrez[!is.na(genes.entrez)]
  x <- enrichPathway(gene=genes.entrez ,pvalueCutoff=0.2, readable=T)
  
  pdf(paste("200703_DSeq2WGCNA/plots/WGCNA/", c, ".reactome.dot.pdf", sep = ""), width = 10, height = 4)
  dotplot(x, orderBy = "x")
  dev.off()
}

for (c in names(table(mergedColors))) {
print(c)
  }
c <- "yellow"
genes <- read.csv(file = paste("200703_DSeq2WGCNA/",c,".txt", sep = "", header=FALSE))
q <- getBM(filters="hgnc_symbol", attributes=c("entrezgene_id", "hgnc_symbol"), values=genes, mart=mart, )
genes.entrez <- q$entrezgene_id
genes.entrez <- genes.entrez[!is.na(genes.entrez)]
x <- enrichPathway(gene=genes.entrez ,pvalueCutoff=0.2, readable=T)

pdf(paste("200703_DSeq2WGCNA/plots/WGCNA/", c, ".reactome.dot.pdf", sep = ""), width = 10, height = 4)
dotplot(x, orderBy = "x")
dev.off()

ego <- enrichGO(gene          = genes.entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dotplot(ego)
cnetplot(ego)

# network analysis for each module
# Select module
module = "yellow";
# Select module probes
genes = colnames(datExpr)
inModule = (moduleColors==module);
modGenes = genes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste(paste("200703_DSeq2WGCNA/TOM", module, ".txt", sep = ""), sep=""),
                            weighted = TRUE,
                            threshold = 0)

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("200703_DSeq2WGCNA/cytoscape/CytoscapeInput-edges-", module, ".txt", sep = ""),
                               nodeFile = paste("200703_DSeq2WGCNA/cytoscape/CytoscapeInput-nodes-", module, ".txt", sep = ""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("200703_DSeq2WGCNA/cytoscape/CytoscapeInput-edges-", "all", ".txt", sep = ""),
                               nodeFile = paste("200703_DSeq2WGCNA/cytoscape/CytoscapeInput-nodes-", "all", ".txt", sep = ""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = genes,
                               nodeAttr = moduleColors)


# Serch module where a gene is included
genes = colnames(datExpr)
moduleColors[genes == "DCLK1"]
moduleColors[genes == "TRPM5"]
moduleColors[genes == "GFI1B"]
moduleColors[genes == "NEFM"]

# PCA for MEs

pca <- plotPCA(vsd, intgroup=c("primary_pathology_histological_type_list"), returnData = TRUE) 
# pca.mean <- dplyr::as_tibble(pca) %>% 
#   group_by(primary_pathology_histological_type_list) %>%
#   summarise(PC1mean = mean(PC1), PC2mean = mean(PC2))
tibble.MEs <- add_rownames(MEs, var = "name")
write.csv(tibble.MEs, "200703_DSeq2WGCNA/ME.sample.csv")
pca <- as_tibble(pca) %>%
  left_join(tibble.MEs, by = "name")

pca.blue <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = MEblue)) +
  ggtitle("MEblue") +
  scale_color_viridis() +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) 
pca.blue
ggsave(filename = "200703_DSeq2WGCNA/plots/WGCNA/PCA.blue.pdf", width = 5, height = 4)

pca.yellow <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = MEyellow)) +
  ggtitle("MEyellow") +
  scale_color_viridis() +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) 
pca.yellow
ggsave(filename = "200703_DSeq2WGCNA/plots/WGCNA/PCA.yellow.pdf", width = 5, height = 4)

pca.turquoise <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = MEturquoise)) +
  ggtitle("MEturquoise") +
  scale_color_viridis() +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) 
pca.turquoise
ggsave(filename = "200703_DSeq2WGCNA/plots/WGCNA/PCA.turquoise.pdf", width = 5, height = 4)

pca.black <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = MEblack)) +
  ggtitle("MEblack") +
  scale_color_viridis() +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) 
pca.black
ggsave(filename = "200703_DSeq2WGCNA/plots/WGCNA/PCA.black.pdf", width = 5, height = 4)

pca.green <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = MEgreen)) +
  ggtitle("MEgreen") +
  scale_color_viridis() +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) 
pca.green
ggsave(filename = "200703_DSeq2WGCNA/plots/WGCNA/PCA.green.pdf", width = 5, height = 4)

pca.red <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = MEred)) +
  ggtitle("MEred") +
  scale_color_viridis() +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) 
pca.red
ggsave(filename = "200703_DSeq2WGCNA/plots/WGCNA/PCA.red.pdf", width = 5, height = 4)

pca.grey <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = MEgrey)) +
  ggtitle("MEgrey") +
  scale_color_viridis() +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) 
pca.grey
ggsave(filename = "200703_DSeq2WGCNA/plots/WGCNA/PCA.grey.pdf", width = 5, height = 4)


pdf(file = "200703_DSeq2WGCNA/plots/WGCNA/PCA.all.pdf", width = 12, height = 10)
gridExtra::grid.arrange(pca.black, pca.turquoise, pca.blue, pca.yellow,
                        pca.green, pca.red, pca.grey, nrow = 3)
dev.off()

# heatmap for selected genes
plot.genes <- gene.names[moduleColors=="yellow"][res[gene.names[moduleColors==module.mg],'log2FoldChange'] > 1]
pdf("200703_DSeq2WGCNA/plots/pheatmap.yellow.lfc1.pdf", width = 16, height = 28)
pheatmap(assay(vsd)[plot.genes,], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")
dev.off()

plot.genes <- gene.names[moduleColors=="yellow"]
pdf("200703_DSeq2WGCNA/plots/pheatmap.yellow.pdf", width = 16, height = 36)
pheatmap(assay(vsd)[plot.genes,], annotation_col=df, color = viridis(10), annotation_colors = ann_colors,
         clustering_distance_rows = "correlation", clustering_method = "ward.D2")
dev.off()

plot.genes <- gene.names[moduleColors=="yellow"]
sample.mg <- row.names(datTraits[datTraits$history_myasthenia_gravis.YES.vs.all == 1,])
plot.col <- colnames(assay(vsd))[order(MEs$MEyellow)]
# plot.col <- c(plot.col[plot.col %in% sample.mg], plot.col[!plot.col %in% sample.mg])
breaksList = seq(-2, 2, by = 0.1)
col.pal <- colorRampPalette(colors = c("blue", "white", "red"))
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- grid::segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.05, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.1, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}
heat <- pheatmap(assay(vsd)[plot.genes,plot.col], 
         annotation_col=df, 
         color = col.pal(length(breaksList)), 
         annotation_colors = ann_colors,
         annotation_names_col = FALSE,
         border_color = NA,
         breaks = breaksList,
         clustering_distance_rows = "correlation", 
         clustering_method = "ward.D2", 
         cluster_cols = FALSE,
         scale = "row",
         show_colnames = FALSE)
keep.genes <- read.csv(file = "200703_DSeq2WGCNA/yellow_selected.txt", header = FALSE, stringsAsFactors = FALSE)$V1
pdf("200703_DSeq2WGCNA/plots/pheatmap.yellow.selectedlabeled.pdf", width = 10, height = 6)
add.flag(heat,
         kept.labels = keep.genes,
         repel.degree = 1)
dev.off()

# krt
plot.genes <- c("KRT23", "KRT18", "KRT36", "KRT17", "KRT7", "KRT1", "KRT4", "KRT38",
"KRT32", "KRT6C", "KRT15", "KRT6A")
pdf("200703_DSeq2WGCNA/plots/pheatmap.krt.pdf", width = 8, height = 4)
pheatmap(assay(vsd)[plot.genes,plot.col], 
                 annotation_col=df, 
                 color = col.pal(length(breaksList)), 
                 annotation_colors = ann_colors,
                 annotation_names_col = FALSE,
                 border_color = NA,
                 breaks = breaksList,
                 clustering_distance_rows = "correlation", 
                 clustering_method = "ward.D2", 
                 cluster_cols = TRUE,
                 scale = "row",
                 show_colnames = FALSE)
dev.off()

tf <- read.csv(file = "data/TF_names_v_1.01.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("TF"))$TF
df.tf <- as.data.frame(row.names(res))
colnames(df.tf) <- c("gene")
df.tf <- transform(df.tf, TF=df.tf$gene %in% tf)
write.csv(df.tf, file = "data/tf.csv")

res.tf <- res[rownames(res) %in% tf, ]
res.tf <- res.tf[order(res.tf$padj),]
res.tf

save.image("/mnt/media32TB/home/yyasumizu/bioinformatics/TCGA_thymoma/200703_DSeq2WGCNA/200703_DESeq2WGCNA.RData")
sessionInfo()
