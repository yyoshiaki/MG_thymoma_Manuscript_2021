library(tidyverse)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(org.Hs.eg.db)

colors <- c("black", 'blue', 'green', 'grey', 'red', 'turquoise', 'yellow')

gc <- list()
for (c in colors) {
  genes <- read.csv(file = paste("200703_DSeq2WGCNA/",c,".txt", sep = ""), header = FALSE)
  ids <- bitr(genes$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  gc[[c]] <- ids$ENTREZID
}

ck <- compareCluster(geneCluster = gc, fun = enrichPathway)
ck <- setReadable(x=ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 

dotplot(ck)

pdf(paste("210617_DSeq2WGCNA_reactome/all.reactome.dot.pdf", sep = ""), width = 16, height = 8)
dotplot(ck)
dev.off()
