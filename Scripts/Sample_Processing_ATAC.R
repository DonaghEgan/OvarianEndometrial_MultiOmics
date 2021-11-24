################################################################################
# Donagh Egan
# POIAZ
# Date: Nov 17th 2021
# 
# Sample: All 
# Description: This script performs the following tasks  
#         1) scATAC-seq processing
#         2) scRNA-seq/scATAC-seq integration
#         3) Peak calling
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
library(ggplot2)
library(Seurat)
library(SingleR)
library(mc2d)
BiocManager::install("SingleR", lib= "/home/degan/R/x86_64-pc-linux-gnu-library/4.1", force = T) 


addArchRGenome("hg38")
addArchRThreads(threads = 1) 
set.seed(123)

# Set up directories and file variables:
################################################################################

SAMPLE.ID <- "ALL"
output.dir <- "/home/degan/Ov_Endo_MultiOmics/Outputs/"
setwd(output.dir)

inputFiles <- list.files(path = "/home/degan/Ov_Endo_MultiOmics/Inputs",
                         pattern = "fragments.tsv.gz$")

inputFiles <- paste("/home/degan/Ov_Endo_MultiOmics/Inputs/", inputFiles, sep="")

sampleNames <- c("3533EL","3571DL","36186L","36639L","366C5L","37EACL","38FE7L","3BAE2L", "3E5CFL", "3CCF1L","3E4D1L")

encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")

# Store patient metadata and colors:

# Make patient sample metadata and color assignments 
sampleColors <- RColorBrewer::brewer.pal(11,"Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors)

# Color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],
                  sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],
                  sampleColors[2],sampleColors[11],sampleColors[1])

sampleAnnot <- data.frame(Sample = c("3533EL","3571DL","36186L","36639L",
                                     "366C5L","37EACL","38FE7L","3BAE2L","3CCF1L","3E4D1L","3E5CFL"),
                          Color = sampleColors,
                          Cancer = c("endometrial","endometrial","endometrial","endometrial","endometrial","endometrial",
                                     "ovarian","ovarian","ovarian","ovarian","ovarian"),
                          Histology = c("endometrioid","endometrioid","endometrioid","endometrioid","endometrioid",
                                        "serous","endometrioid","serous","carcinosarcoma","GIST","serous"),
                          BMI = c(39.89,30.5,38.55,55.29,49.44,29.94,34.8,22.13,23.72,33.96,22.37),
                          Age = c(70,70,70,49,62,74,76,61,69,59,59),
                          Race = c("AA","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","AS"),
                          Stage = c("IA","IA","IA","IA","IA","IIIA","IA","IIB","IVB","IV","IIIC"),
                          Site = c("Endometrium","Endometrium","Endometrium","Endometrium","Endometrium","Ovary","Ovary","Ovary","Ovary","Ovary","Ovary"),
                          Type = c("Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Ovarian","Ovarian","Ovarian","Gastric","Ovarian"))

# Create Arrow and ArchR project
################################################################################
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = T,
  addGeneScoreMat = T,
)

ArrowFiles <- list.files(pattern=".arrow", path = "/home/degan/Ov_Endo_MultiOmics/Inputs/")

setwd("/home/degan/Ov_Endo_MultiOmics/Inputs/")
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP",useMatrix = "TileMatrix",nTrials=5,LSIMethod = 1,scaleDims = F,
  corCutOff = 0.75,UMAPParams = list(n_neighbors =30, min_dist = 0.3, metric = "cosine", verbose =FALSE),
  dimsToUse = 1:50
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/home/degan/Ov_Endo_MultiOmics/Outputs/",
  copyArrows = T #This is recommened so that you maintain an unaltered copy for later usage.
)

# Filter out outlier low quality cells and doublets
################################################################################
# GMM for fragments per cell

setwd("/home/degan/Ov_Endo_MultiOmics/Outputs/")

for (i in sampleNames){
  proj.i <- proj[proj$Sample == i]
  
  # GMM for fragments per cell
  depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
  proj.i$depth.cluster <- depth.clust$classification
  proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$depth.cluster),
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))
    ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
  
  # GMM for TSS per cell
  TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
  proj.i$TSS.cluster <- TSS.clust$classification
  proj.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$TSS.cluster),
    discrete = T,
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))
    ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
  
  
  df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
  saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
  
  df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
  saveRDS(df.depth,paste0("df_depth_",i,".rds"))
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n",i))
    ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = proj.i$DoubletEnrichment,
    discrete = F,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("Doublet Enrichment:\n",i))
    ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
  
}


for (i in sampleNames[2]){
  proj.i <- proj[proj$Sample == i]
  
  # GMM for fragments per cell
  depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
  proj.i$depth.cluster <- depth.clust$classification
  proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$depth.cluster),
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))
    ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
  
  # Manually set TSS threshold
  #TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
  proj.i$TSS.cluster <- ifelse(log10(proj.i$TSSEnrichment+1) >= 0.80,"2","1")
  proj.i$TSS.cluster.uncertainty <- rep(NA,nrow(proj.i@cellColData))
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$TSS.cluster),
    discrete = T,
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))
    ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
  
  
  df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
  #df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
  saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
  
  df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
  saveRDS(df.depth,paste0("df_depth_",i,".rds"))
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n",i))
    ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = proj.i$DoubletEnrichment,
    discrete = F,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("Doublet Enrichment:\n",i))
    ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
  
}

for (i in sampleNames[7]){
  proj.i <- proj[proj$Sample == i]
  
  # GMM for fragments per cell
  depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
  proj.i$depth.cluster <- depth.clust$classification
  proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$depth.cluster),
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))
    ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
  
  # Manually set TSS threshold
  #TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
  proj.i$TSS.cluster <- ifelse(log10(proj.i$TSSEnrichment+1) >= 0.80,"2","1")
  proj.i$TSS.cluster.uncertainty <- rep(NA,nrow(proj.i@cellColData))
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$TSS.cluster),
    discrete = T,
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))
    ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
  
  
  df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
  #df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
  saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
  
  df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
  saveRDS(df.depth,paste0("df_depth_",i,".rds"))
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n",i))
    ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = proj.i$DoubletEnrichment,
    discrete = F,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("Doublet Enrichment:\n",i))
    ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
  
}

################################################################################
dev.off()

# Filter out low quality cells, and remove doublets
################################################################################

list.depth <- list.files(pattern = "^df_depth")

df.depth <-  data.frame(cellNames=character(),
                        cluster=character(),
                        cluster.uncertainty=character(),
                        nFrags = character())
for (i in list.depth){
  df <- readRDS(i)
  colnames(df) <- c("cellNames","cluster","cluster.uncertainty","nFrags")
  df.depth <- rbind(df.depth,df)
}

list.TSS <- list.files(pattern = "^df_TSS")

df.TSS <-  data.frame(cellNames=character(),
                      cluster=character(),
                      cluster.uncertainty=character(),
                      TSSEnrichment = character())
for (i in list.TSS){
  df <- readRDS(i)
  colnames(df) <- c("cellNames","cluster","cluster.uncertainty","TSSEnrichment")
  df.TSS <- rbind(df.TSS,df)
}


colnames(df.TSS) <- c("cellNames","TSS.cluster","TSS.cluster.uncertainty","TSSEnrichment")
colnames(df.depth) <- c("cellNames","depth.cluster","depth.cluster.uncertainty","nFrags")

cellsPass <- intersect(df.TSS$cellNames,df.depth$cellNames)

cellsFail <-  proj$cellNames[!(proj$cellNames %in% cellsPass)]

# Screen for high quality barcodes (remove non cellular barcodes)
proj.filter <- proj[proj$cellNames %in% cellsPass]


proj <- filterDoublets(proj.filter,filterRatio = 1,cutEnrich = 1,cutScore = -Inf)


plotFragmentSizes(proj)+ggtitle("Fragment Size Histogram")
ggsave("Frags_hist.pdf",width = 6,height = 4)
plotTSSEnrichment(proj)+ggtitle("TSS Enrichment")
ggsave("TSS.pdf",width = 6,height = 4)

################################################################################


# Perform LSI reduction and clustering with ATAC data only
################################################################################

# Add LSI dimreduc
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  LSIMethod = 2,
  scaleDims = T,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  UMAPParams = list(n_neighbors =30,
                    min_dist = 0.3,
                    metric = "cosine",
                    verbose =FALSE),
  varFeatures = 25000,
  dimsToUse = 1:50,
  binarize = T,
  corCutOff = 0.75,
  force = T,
  seed=6
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ATAC_clusters",
  resolution = 0.7,
  dimsToUse = 1:50,force = T
)

# Add UMAP based on LSI dims
proj <- addUMAP(proj,nNeighbors = 30,minDist = 0.3,dimsToUse = 1:50,metric = "cosine",force = T,reducedDims="IterativeLSI")
saveRDS(proj,"proj_LSI_AND_UMAP.rds")

################################################################################

# Estimate gene activity in ATAC data and perform cell type annotation:

# Add Gene activity matrix using ArchR model
proj <- addGeneScoreMatrix(proj,matrixName = "ArchRGeneScore",force = T)


getAvailableMatrices(proj)
saveRDS(proj,"proj_LSI_AND_GeneScores.rds")


################################################################################
# Begin Downstream Analysis
#      1)  Plotting RNA/ATAC by sample, by cluster, by predicted label
#      2)  Marker Gene (RNA/ATAC) intersection
#      3)  Peak2GeneLinks/Coaccessiblity
######################################################

# PART 1: Plotting
######################################################################################

# Make embedding highlighting by 1) Predicted group ArchR 2) Predicted group Signac
# 3) Sample 4) ATAC-only clusters
atac <- plotEmbedding(proj,colorBy = "cellColData",name = "ATAC_clusters")
atac.emb.cluster <- as.data.frame(atac$data)
atac.emb.cluster$sample <- atac.emb.cluster$color
atac.emb.cluster$sample <- sub("-", ":", atac.emb.cluster$sample )
atac.emb.cluster$sample  <- gsub(".*:", "", atac.emb.cluster$sample )
head(atac.emb.cluster)

atac <- plotEmbedding(proj,colorBy = "cellColData",name = "Sample")
atac.emb.sample <- as.data.frame(atac$data)
atac.emb.sample$sample <- atac.emb.sample$color
atac.emb.sample$sample <- sub("-", ":", atac.emb.sample$sample )
atac.emb.sample$sample  <- gsub(".*:", "", atac.emb.sample$sample )
head(atac.emb.sample)

atac.emb.all <- cbind(atac.emb.sample[,c(1:2,4)],atac.emb.cluster[,4])
names(atac.emb.all)[4] <- "cluster"
atac.emb.all$plain <- "Plain"

colnames(atac.emb.all) <- c("UMAP1","UMAP2","Sample", "ATAC_clusters", "Blank")

ggplot(atac.emb.all,aes_string(x = "UMAP1",y="UMAP2",color = "Sample"))+
  geom_point(size = .1)+
  theme_classic()+
  ggtitle(paste0("scATAC-seq: Sample"))+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)
#guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)