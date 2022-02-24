################################################################################
# Donagh Egan
# POIAZ 
# Date: 24 December 2021
# 

# Description: This script combines and processes the individual counts data
#              for patients 1-6 
#         
################################################################################


## library and internal functions ####
################################################################################
source("/home/degan/Ov_Endo_MultiOmics/Scripts/rowr.R")
source("/home/degan/Ov_Endo_MultiOmics/Scripts/stacked_violin.R")

library(scater)
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(scater)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tibble)
library(DoubletFinder)
library(Signac)
library(ggplot2)
library(stringr)
library(SingleR)
library(psych)

## loading gene signatures ####
################################################################################

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

ESTIMATE.signatures <- "./ESTIMATE_signatures.csv"

SAMPLE.ID = "endo_ovar_All"

## Part 1: Merge individually processed scRNA-seq datasets ####
################################################################################

# list processed scRNA-seq Seurat objects made in the previous R scripts
datasets <- list.files(pattern = "*_processed.rds")

endo_3533EL <- readRDS(datasets[1])
endo_3533EL$SingleR <- endo_3533EL$SingleR.endo
endo_3533EL$Sample <- "3533EL"

endo_ovar_37EACL <- readRDS(datasets[2])
endo_ovar_37EACL$SingleR <- endo_ovar_37EACL$SingleR.endo
endo_ovar_37EACL$Sample <- "37EACL"

# Merge Seurat objects
rna <- merge(x = endo_3533EL,
             y = list(endo_ovar_37EACL))

rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna)
rna <- ScaleData(rna,features = all.genes)
rna <- RunPCA(rna)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")

## Score cells for immune/stromal/fibroblast/endothelial signatures ####
################################################################################


stromal <- ESTIMATE.signatures$X2[3:143]
immune <- ESTIMATE.signatures$X3[3:143]
fibroblast <- panglaodb$Fibroblasts
endothelial <- panglaodb$`Endothelial cells`
epithelial <- panglaodb$`Epithelial cells`
smooth <- panglaodb$`Smooth muscle cells`
plasma <- panglaodb$`Plasma cells`

feats <- list(stromal,immune,fibroblast,endothelial,epithelial,smooth,plasma)

rna <- AddModuleScore(rna,features = feats,name = c("stromal.","immune.","fibroblast.","endothelial.",
                                                    "epithelial.","smooth.","plasma."),search = T)


## Add PC1 to metadata and check correlations ####
################################################################################
rna@meta.data$PC1 <- rna@reductions$pca@cell.embeddings[,1]

count_cor_PC1 <- cor(rna$PC1,rna$nCount_RNA,method = "spearman")

stromal.cor <- cor(rna$PC1,rna$stromal.1,method = "spearman")
immune.cor <- cor(rna$PC1,rna$immune.2,method = "spearman")
fibroblast.cor <- cor(rna$PC1,rna$fibroblast.3,method = "spearman")
endothelial.cor <- cor(rna$PC1,rna$endothelial.4,method = "spearman")
epithelial.cor <- cor(rna$PC1,rna$epithelial.5,method = "spearman")
smooth.cor <- cor(rna$PC1,rna$smooth.6,method = "spearman")
plasma.cor <- cor(rna$PC1,rna$plasma.7,method = "spearman")

###########################################################
## Part 2: Reprocessing and CNV-read depth check #### 
###########################################################

# If PC1 is correlated with read depth, check to see if biological variation is correlated to PC1
if (round(abs(count_cor_PC1),2) > 0.5){
  
  if( round(abs(stromal.cor),2) >= 0.5 |
      round(abs(immune.cor),2) >= 0.5 |
      round(abs(fibroblast.cor),2) >= 0.5 |
      round(abs(endothelial.cor),2) >= 0.5 |
      round(abs(epithelial.cor),2) >= 0.5 |
      round(abs(smooth.cor),2) >= 0.5 |
      round(abs(plasma.cor),2) >= 0.5){
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    
    rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
    rna$cell.barcode <- rownames(rna@meta.data)
    rna$CNV.Pos <- ifelse(rna$Total_CNVs > 0,TRUE,FALSE)
    
    boxplot.cnv <- ggplot(rna@meta.data,aes(x= RNA_snn_res.0.7,y=PC1.loading,color = as.factor(CNV.Pos)))+geom_boxplot()
    boxplot.cnv
    ggsave("CNV_PC1_boxplot.pdf", width = 10, height = 4)
    
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

###########################################################
## Part 3: Cell type annotation of clusters ####
###########################################################

Idents(rna) <- "RNA_snn_res.0.7"

DimPlot(rna,group.by = "Sample",label = T)
ggsave("Patient_UMAP.pdf",width = 6,height = 4)
DimPlot(rna,group.by = "RNA_snn_res.0.7",label = T)
ggsave("Cluster_UMAP.pdf",width = 6,height = 4)
DimPlot(rna,group.by = "SingleR",label = T)
ggsave("SingleR_UMAP.pdf",width = 6,height = 4)

# Verify SingleR annotations and check for mast cells:
rna <- AddModuleScore(rna,features = list(panglaodb$`B cells`,
                                          panglaodb$`Plasma cells`,
                                          panglaodb$`Mast cells`,
                                          panglaodb$Macrophages,
                                          panglaodb$`Dendritic cells`,
                                          panglaodb$`T cells`,
                                          panglaodb$`NK cells`,
                                          panglaodb$`Endothelial cells`,
                                          panglaodb$Fibroblasts,
                                          panglaodb$`Epithelial cells`,
                                          panglaodb$`Smooth muscle cells`,
                                          c("TPSB2","TPSAB1","KIT"), #Three gene Mast signature
                                          c("MUC16", "WFDC2")), #cancer genes 
                      name = c("B.","Plasma.","Mast.","Macrophage.","DC.",
                               "T.","NK.","Endothelial.","Fibroblast.","Epithelial.","Smooth_muscle.","Mast_3_gene.", "Cancer Cells"),search = T)

## Visualise gene signatures 
StackedVlnPlot(rna,features = c("B.1","Plasma.2","Mast.3","Macrophage.4",
                                "DC.5","T.6","NK.7","Endothelial.8","Fibroblast.9",
                                "Epithelial.10","Smooth_muscle.11"))
ggsave("Panglaodb_Signatures_Violin.pdf",width = 8,height = 16)

StackedVlnPlot(rna,features = c("TPSB2","TPSAB1","KIT"))
ggsave("Mast_Signatures_Violin.pdf",width = 8,height = 4)

StackedVlnPlot(rna,c("MUC16","WFDC2"))
ggsave("Cancer_Signatures_Violin.pdf",width = 8,height = 3)

## Assess mast cell, bcell and cancer cell enrichment to potentially rename clusters ####
#########################################################################################

## Assess Cancer cell enrichment 
vln.df <- VlnPlot(rna,features = "Cancer.Cells13")
data.cancer <- describeBy(vln.df$data$Cancer.Cells13, vln.df$data$ident, mat = TRUE)
data.cancer <- dplyr::filter(data.cancer,median > 0.1)

## Assess Mast cell enrichment 
vln.df <- VlnPlot(rna,features = "Mast_3_gene.12")
data.mast <- describeBy(vln.df$data$Mast_3_gene.12, vln.df$data$ident, mat = TRUE)
data.mast <- dplyr::filter(data.mast,median > 0.225)

## Assess B cell enrichment 
vln.df <- VlnPlot(rna,features = "B.1")
data.B <- describeBy(vln.df$data$B.1, vln.df$data$ident, mat = TRUE)
data.B <- dplyr::filter(data.B,median > 0.225)

# Annotate mast/b/cancer cells
rna$mast.cell <- ifelse(rna$RNA_snn_res.0.7 %in% as.character(data.mast$group1),TRUE,FALSE)
rna$B.cell <- ifelse(rna$RNA_snn_res.0.7 %in% as.character(data.B$group1),TRUE,FALSE)
rna$Cancer.cell <- ifelse(rna$RNA_snn_res.0.7 %in% as.character(data.cancer$group1),TRUE,FALSE)

cells <- rna@meta.data %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(SingleR) 

cluster.ids <- rep("fill",length(levels(Idents(rna))))
names(cluster.ids) <- levels(Idents(rna))

for ( i in factor(cells$RNA_snn_res.0.7)){
  library(tidyr)
  cells.sub <- cells %>% dplyr::filter(RNA_snn_res.0.7 ==i) %>% arrange(desc(n))
  cluster.ids[[i]] <- cells.sub$SingleR[1]
  
}

# Rename cluster if median enrichment score is greater than 0.1  
if(nrow(data.mast) > 0){
  for (i in 1:nrow(data.mast)){
    cluster.ids[[data.mast$group1[i]]] <- "Mast cell" #Marker Mast cell cluster 
  }
}else{
  cluster.ids <- cluster.ids
}

# Rename cluster if median enrichment score is greater than 0.1  
if (nrow(data.B) >0 ){
  for (i in 1:nrow(data.B)){
    cluster.ids[[data.B$group1[i]]] <- "B cell" #Marker B cell cluster 
  }
}else{
  cluster.ids <- cluster.ids
}


# Rename cluster if median enrichment score is greater than 0.1  
if (nrow(data.cancer) >0 ){
  for (i in 1:nrow(data.cancer)){
    cluster.ids[[data.cancer$group1[i]]] <- "Cancer cell" #Marker Cancer cell cluster 
  }
}else{
  cluster.ids <- cluster.ids
}

cluster.ids <- as.data.frame(cluster.ids)

levels(Idents(rna)) <- cluster.ids$cluster.ids

rna$cell.type <- Idents(rna)

rna$cell.type <- paste0(rna$RNA_snn_res.0.7,"-",rna$cell.type)

Idents(rna) <- rna$cell.type

DimPlot(rna,group.by = "cell.type",label = F)
ggsave("cell_type_UMAP.pdf",width = 8,height = 5)



























