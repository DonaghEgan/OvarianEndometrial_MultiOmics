## Reading in and analyzing data from Izar et al 2020
## https://www.nature.com/articles/s41591-020-0926-0#Abs1



## Reading in ov counts ####
## 
################################################################################


## ss2 = SmartSeq2 - used to sequence epithelial cells
ov_counts_ss2 <- read.table(file = "/home/degan/Ov_Endo_MultiOmics/Inputs/ref_datasets/GSE146026_Izar_HGSOC_ascites_SS2_log.tsv", header = TRUE)

ov_counts <- data.frame(ov_counts_ss2[c(6:length(rownames(ov_counts_ss2))),])
names(ov_counts)[1] <- "Gene"
rownames(ov_counts) <- NULL
ov_counts <- column_to_rownames(ov_counts, "Gene")

col_data <- data.frame(ov_counts_ss2[c(1:5),])
col_data <- column_to_rownames(col_data, "Cell_ID")
col_data <- data.frame(t(col_data))
col_data$celltype <- ifelse(col_data$clst %in% (1:5), "Malignant", "Fibroblast")

sce <- SingleCellExperiment(list(counts=ov_counts),
                            colData=DataFrame(label=col_data),
                            metadata=list(study="GSE146026")
)
