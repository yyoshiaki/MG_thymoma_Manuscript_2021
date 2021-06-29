library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)
library(TCGAmutations)
library(dplyr)
library(DT)

tcga_load(study = "THYM")
getSampleSummary(x = tcga_thym_mc3)

getGeneSummary(x = tcga_thym_mc3)

getClinicalData(x = tcga_thym_mc3)[1:10, 1:10]

frame()
plotmafSummary(maf = tcga_thym_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plot.new()  
oncoplot(maf = tcga_thym_mc3) #, top = 5, fontSize = 12)

plot.new()  
oncoplot(maf = tcga_thym_mc3, clinicalFeatures = 'primary_pathology_history_myasthenia_gravis', sortByAnnotation = TRUE)

thym.sig = oncodrive(maf = tcga_thym_mc3, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

