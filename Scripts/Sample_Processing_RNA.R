################################################################################
# Donagh Egan
# POIAZ
# Date: Nov 17th 2021
# 
# Sample: All 
# Description: This script performs the following tasks  
#         1) scRNA-seq processing
################################################################################

# Library 
################################################################################
library(mclust)
library(RColorBrewer)
library(chromVAR)
library(ArchR)
library(Matrix)
library(dplyr)
library(tibble)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(ensembldb)
library(readxl)
library(ggplot2)
library(stringr)
library(scater)

BiocManager::install("scater", lib= "/home/degan/R/x86_64-pc-linux-gnu-library/4.1", force = T) 


set.seed(123)

## Set up directories and file variables:
################################################################################
setwd("/home/degan/Ov_Endo_MultiOmics/Inputs/")

sampleNames <- c("3533EL","3571DL","36186L","36639L","366C5L","37EACL","38FE7L","3BAE2L", "3CCF1L","3E4D1L", "3E5CFL")


inputRNA <- list.files(path = "/home/degan/Ov_Endo_MultiOmics/Inputs",
                         pattern = "mtx.gz$")

barcodefiles <- list.files(path = "/home/degan/Ov_Endo_MultiOmics/Inputs",
                           pattern = "barcodes")

featurefiles <- list.files(path = "/home/degan/Ov_Endo_MultiOmics/Inputs",
                           pattern = "features")

## PanglaoDB 
tsv=gzfile("./PanglaoDB_markers_27_Mar_2020.tsv.gz")  
panglaodb <- read.csv(tsv,header=T,sep = "\t") 
panglaodb <- dplyr::filter(panglaodb,species == "Hs" | species == "Mm Hs")# Human subset 
panglaodb <- dplyr::filter(panglaodb,organ == "Connective tissue" |
                             organ == "Epithelium" |
                             organ == "Immune system" |
                             organ == "Reproductive"|
                             organ == "Vasculature" |
                             organ == "Smooth muscle"
)
panglaodb <- split(as.character(panglaodb$official.gene.symbol), panglaodb$cell.type)

## ESTIMATE signature - stromal and immune 
ESTIMATE.signatures <- read_excel("41467_2013_BFncomms3612_MOESM488_ESM.xlsx")
ESTIMATE.signatures <- data.frame(t(ESTIMATE.signatures))

## Load in RNAseq counts and merge into a seurat object 
seurat_objects <-list()
for (i in c(1:11)) {
  RNAcounts <- ReadMtx(mtx = inputRNA[i], cells = barcodefiles[i], features = featurefiles[i])
  RNAcounts_seurat <- CreateSeuratObject(counts = RNAcounts, min.cells = 3)
  seurat_objects[i] <- RNAcounts_seurat
}
names(seurat_objects) <- sampleNames

seurat_combined <- merge(x=seurat_objects[["3533EL"]],
                         y= c(seurat_objects[["3571DL"]],
                              seurat_objects[["36186L"]],
                              seurat_objects[["36639L"]],
                              seurat_objects[["366C5L"]],
                              seurat_objects[["37EACL"]],
                              seurat_objects[["38FE7L"]],
                              seurat_objects[["3BAE2L"]],
                              seurat_objects[["3CCF1L"]],
                              seurat_objects[["3E4D1L"]],
                              seurat_objects[["3E5CFL"]]))

PreQCNumCells <- length(colnames(seurat_combined))

## Outlier detection: nCounts, nfeatures, and percent mitochondrial counts
## Outliers are >2 MADs (2 SDs from the median)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")

seurat_combined@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(seurat_combined@meta.data$nCount_RNA),
                                                              log = F,type = "lower",nmads = 2)
seurat_combined@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(seurat_combined@meta.data$nFeature_RNA),
                                                                  log = F,type = "lower",nmads = 2)
seurat_combined@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(seurat_combined@meta.data$percent.mt),
                                                   log = F,type = "higher",nmads = 2)

seurat_combined <- subset(seurat_combined, subset = nCount_RNA_outlier_2mad == "FALSE" &
                nFeature_RNA_outlier_2mad == 'FALSE' & 
                percent_mt_outlier_2mad == "FALSE")

PostQCNumCells <- length(colnames(seurat_combined))


seurat_combined <- NormalizeData(seurat_combined, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_combined)
seurat_combined <- ScaleData(seurat_combined,features = all.genes)
seurat_combined <- RunPCA(seurat_combined)

# Score cells for immune/stromal/fibroblast/endothelial signatures
############################################################################################


stromal <- ESTIMATE.signatures$X2[3:143]
immune <- ESTIMATE.signatures$X3[3:143]
fibroblast <- panglaodb$Fibroblasts
endothelial <- panglaodb$`Endothelial cells`
epithelial <- panglaodb$`Epithelial cells`
smooth <- panglaodb$`Smooth muscle cells`
plasma <- panglaodb$`Plasma cells`

feats <- list(stromal,immune,fibroblast,endothelial,epithelial,smooth,plasma)

seurat_combined <- AddModuleScore(seurat_combined,features = feats,name = c("stromal.","immune.","fibroblast.","endothelial.",
                                                    "epithelial.","smooth.","plasma."),search = T)

#########################################################################

#######################################################################

## Add PC1 to metadata
seurat_combined@meta.data$PC1 <- seurat_combined@reductions$pca@cell.embeddings[,1]

count_cor_PC1 <- cor(seurat_combined$PC1, seurat_combined$nCount_RNA, method = "spearman")

## checking correlation between cell signature scores and PC1 
stromal.cor <- cor(seurat_combined$PC1,seurat_combined$stromal.1,method = "spearman")
immune.cor <- cor(seurat_combined$PC1,seurat_combined$immune.2,method = "spearman")
fibroblast.cor <- cor(seurat_combined$PC1,seurat_combined$fibroblast.3,method = "spearman")
endothelial.cor <- cor(seurat_combined$PC1,seurat_combined$endothelial.4,method = "spearman")
epithelial.cor <- cor(seurat_combined$PC1,seurat_combined$epithelial.5,method = "spearman")
smooth.cor <- cor(seurat_combined$PC1,seurat_combined$smooth.6,method = "spearman")
plasma.cor <- cor(seurat_combined$PC1,seurat_combined$plasma.7,method = "spearman")

###########################################################
# Part 2: Reprocessing and CNV-read depth check
###########################################################

# If PC1 is correalted with read depth, check to see if biological variation is corralted to PC1
if (round(abs(count_cor_PC1),2) > 0.5){
  
  if( round(abs(stromal.cor),2) >= 0.5 |
      round(abs(immune.cor),2) >= 0.5 |
      round(abs(fibroblast.cor),2) >= 0.5 |
      round(abs(endothelial.cor),2) >= 0.5 |
      round(abs(epithelial.cor),2) >= 0.5 |
      round(abs(smooth.cor),2) >= 0.5 |
      round(abs(plasma.cor),2) >= 0.5){
    seurat_combined <- FindNeighbors(seurat_combined,dims = 1:50)
    seurat_combined <- FindClusters(seurat_combined,resolution = 0.7)
    seurat_combined <- RunUMAP(seurat_combined,dims = 1:50)
    Idents(seurat_combined) <- "RNA_snn_res.0.7"

    rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
    rna$cell.barcode <- rownames(rna@meta.data)
    rna$CNV.Pos <- ifelse(rna$Total_CNVs > 0,TRUE,FALSE)
    
    boxplot.cnv <- ggplot(rna@meta.data,aes(x= RNA_snn_res.0.7,y=PC1.loading,color = as.factor(CNV.Pos)))+geom_boxplot()
    boxplot.cnv+ggsave("CNV_PC1_boxplot.png")
    
    data <- describeBy(boxplot.cnv$data$PC1.loading, boxplot.cnv$data$RNA_snn_res.0.7, mat = TRUE)
    
    wilcox <- wilcox.test(data = rna@meta.data,PC1.loading~CNV.Pos)
    
    if (wilcox$p.value < 0.05){
      rna <- rna
      saveRDS(rna,"./rna_PassedPC1Checks.rds")
    }else{
      all.genes <- rownames(rna)
      rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
      rna <- RunPCA(rna)
      rna <- FindNeighbors(rna,dims = 1:50)
      rna <- FindClusters(rna,resolution = 0.7)
      rna <- RunUMAP(rna,dims = 1:50)
      Idents(rna) <- "RNA_snn_res.0.7"
      saveRDS(rna,"./rna_FailedCNVTest.rds")
    }
    
  }else{
    all.genes <- rownames(rna)
    rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
    rna <- RunPCA(rna)
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    saveRDS(rna,"./rna_FailedCorTest.rds")
  }
}else{
  rna <- FindNeighbors(rna,dims = 1:50)
  rna <- FindClusters(rna,resolution = 0.7)
  rna <- RunUMAP(rna,dims = 1:50)
  Idents(rna) <- "RNA_snn_res.0.7"
  saveRDS(rna,"./rna_SkipChecks.rds")
}





