options(scipen = 999)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("recount3")
#BiocManager::install("recount")
#BiocManager::install("sva")

library("recount3")
library("recount")
library("data.table")
library(Biobase)
library(readr)
library(rtracklayer)
library(biomaRt)
library(edgeR) ## For TMM
library(sva) ## For batch correction
library(dplyr)

GTEx = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/GTEx_RSE.RData"))
TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))

# GTEx & TCGA raw counts
assays(GTEx_lung)$counts <- recount3::transform_counts(GTEx_lung)
assays(TCGA_lung)$counts <- recount3::transform_counts(TCGA_lung)

##### Convert raw counts into TPM
assays(GTEx_lung)$TPM = recount::getTPM(GTEx_lung)
assays(TCGA_lung)$TPM = recount::getTPM(TCGA_lung)

GTEx_tpm = assays(GTEx_lung)$TPM
TCGA_tpm = assays(TCGA_lung)$TPM

# Remove GTEx samples not present in portal
gtex_portal_id <- colnames(read.csv("/home/ubuntu/lung_expression_v8.txt", sep=""))
gtex_portal_id = unlist(lapply(strsplit(gtex_portal_id, split = ".", fixed = T), function(x){paste(paste(x, collapse = "-"), 1, sep = ".")}))

GTEx_tpm = GTEx_tpm[,which(colnames(GTEx_tpm) %in% gtex_portal_id)]

# Extract primary tumor samples from TCGA
TCGA_tpm = TCGA_tpm[,which(TCGA_lung@colData$tcga.cgc_sample_sample_type == "Primary Tumor")]

# Keep the sample with maximum depth set and remove other duplicate samples
remove_duplicate = function(expression_matrix, unique_ids){
  removed_samples = c()
  for (i in 1:length(unique(unique_ids))){
    id = unique(unique_ids)[i]
    id_loc = which(unique_ids == id)
    if (length(id_loc) > 1){
      subdata = expression_matrix[,id_loc]
      seq_depth = colSums(subdata)
      removed_samples = c(removed_samples, colnames(subdata)[-which.max(seq_depth)])
    }
  }
  expression_matrix = expression_matrix[,-which(colnames(expression_matrix) %in% removed_samples)]
  return(expression_matrix)
}

# get donor & tissue ids of GTEx:
GTEx_ID = unlist(lapply(strsplit(colnames(GTEx_tpm), split = "-"), function(x){paste(x[1:2], collapse = "-")}))
# get donor ids for TCGA
TCGA_ID = sapply(strsplit(TCGA_lung@colData$tcga.tcga_barcode[which(TCGA_lung@colData$tcga.cgc_sample_sample_type == "Primary Tumor")], split="-"),
                 function(x){paste(x[1], x[2], x[3], x[4], sep = "-")})
# remove vial (last letter) to identify unique donors
TCGA_ID = substr(TCGA_ID, 1,nchar(TCGA_ID)-1)

# Samples with duplicates
which(table(GTEx_ID) > 1)
which(table(TCGA_ID) > 1)

#GTEx_tpm = remove_duplicate(GTEx_tpm, GTEx_ID)
TCGA_tpm = remove_duplicate(TCGA_tpm, TCGA_ID)

##### Remove genes with counts <1 TPM in at least 10% samples
expression_all = cbind(GTEx_tpm,TCGA_tpm)
expression_cutoff = 1 # 1 TPM
percent_cutoff = 0.1
minSamples = percent_cutoff*ncol(expression_all) # at least 5% of samples
keep = rowSums(expression_all > expression_cutoff) >= minSamples
table(keep)
expression_filtered = expression_all[keep,]
paste0(nrow(expression_all), ":Total number of genes overlapped between TCGA and GTEx")
paste0(nrow(expression_filtered), ":Number of genes after filtering ", expression_cutoff," cpm in ", minSamples, " samples")

# Add filtered TPM to data
GTEx_TPM = expression_filtered[,1:ncol(GTEx_tpm)]
TCGA_TPM = expression_filtered[,-(1:ncol(GTEx_tpm))]

# Transform to log2(1+TPM)
GTEx_TPM = log2(1+GTEx_TPM)
TCGA_TPM = log2(1+TCGA_TPM)

#######################################
# Get phenotypic data
GTEx_phenotypes = GTEx_lung@colData[match(colnames(GTEx_TPM),rownames(GTEx_lung@colData)),]
TCGA_phenotypes = TCGA_lung@colData[match(colnames(TCGA_TPM),rownames(TCGA_lung@colData)),]

# Get GTEx phenotypes
# Remove cancer samples from GTEx
gtex.v8_phenotypes = data.frame(fread("/home/ubuntu/GTEx.v8.phenotypes.txt", skip = 10))
GTEx_phenotypes = cbind(GTEx_phenotypes, gtex.v8_phenotypes[match(GTEx_phenotypes$gtex.subjid, gtex.v8_phenotypes$SUBJID),])

# PCA of Y genes to identify mismatch between gender and sex
# Get chromosome location
fdata.GTEx = as.data.frame(GTEx_lung@rowRanges)
fdata.TCGA = as.data.frame(TCGA_lung@rowRanges)

### Get GenecIDs for Y genes & XIST
Y_genes.GTEx = fdata.GTEx$gene_id[which(fdata.GTEx$seqnames == "chrY")]
Y_genes.TCGA = fdata.TCGA$gene_id[which(fdata.TCGA$seqnames == "chrY")]

# Get Y gene expression
Y_gtex = GTEx_TPM[which(rownames(GTEx_TPM) %in% Y_genes.GTEx),]
Y_tcga = TCGA_TPM[which(rownames(TCGA_TPM) %in% Y_genes.TCGA),]

# PCA of Y genes
library(ggplot2)

pca.GTEx <- prcomp(t(Y_gtex))
U <- pca.GTEx$x
U.GTEx <- U
GTEx_sex = GTEx_phenotypes$gtex.sex
GTEx_sex[which(GTEx_sex == 1)] = "MALE"
GTEx_sex[which(GTEx_sex == 2)] = "FEMALE"
GTEx_phenotypes$SEX = GTEx_sex
cols = factor(GTEx_sex)
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], main = "GTEx: PCA of Y genes", colour = cols, xlab = "PC1", ylab = "PC2") + coord_equal()

pca.TCGA <- prcomp(t(Y_tcga))
U <- pca.TCGA$x
U.TCGA <- U
cols = factor(TCGA_phenotypes$tcga.cgc_case_gender)
TCGA_sex = cols
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], main = "TCGA: PCA of Y genes", colour = cols, xlab = "PC1", ylab = "PC2") + coord_equal()

# Find sex-misannotated samples
# GTEx_misannotated = colnames(GTEx_TPM)[which(U.GTEx[,1] < 0 & GTEx_sex == "FEMALE")]
TCGA_misannotated = colnames(TCGA_TPM)[intersect(which(U.TCGA[,1] > 0), which(TCGA_sex == "FEMALE"))]
# exclude 1 sample with gender = "NA" in TCGA
TCGA_misannotated  = c(TCGA_misannotated, rownames(TCGA_phenotypes)[is.na(TCGA_phenotypes$tcga.cgc_case_gender)])

# Exclude samples with sex misannotation
#GTEx_TPM = GTEx_TPM[,-which(colnames(GTEx_TPM) %in% GTEx_misannotated)]
TCGA_TPM = TCGA_TPM[,-which(colnames(TCGA_TPM) %in% TCGA_misannotated)]

# Get phenotypic data
GTEx_phenotypes = GTEx_phenotypes[match(colnames(GTEx_TPM),rownames(GTEx_phenotypes)),]
TCGA_phenotypes = TCGA_phenotypes[match(colnames(TCGA_TPM),rownames(TCGA_phenotypes)),]

write.table(GTEx_TPM, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_TPM.txt", row.names = T, col.names = T)
write.table(TCGA_TPM, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_TPM.txt", row.names = T, col.names = T)
write.table(GTEx_phenotypes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_phenotypes.txt", row.names = T, col.names = T)
write.table(TCGA_phenotypes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_phenotypes.txt", row.names = T, col.names = T)

##############
# Separate Male and Female Expression
male.GTEx = colnames(GTEx_TPM)[which(GTEx_phenotypes$SEX == "MALE")]
female.GTEx = colnames(GTEx_TPM)[which(GTEx_phenotypes$SEX == "FEMALE")]
male.TCGA = colnames(TCGA_TPM)[which(TCGA_phenotypes$tcga.cgc_case_gender == "MALE")]
female.TCGA = colnames(TCGA_TPM)[which(TCGA_phenotypes$tcga.cgc_case_gender == "FEMALE")]

GTEx_TPM.male = GTEx_TPM[,which(colnames(GTEx_TPM) %in% male.GTEx)]
GTEx_TPM.female = GTEx_TPM[,which(colnames(GTEx_TPM) %in% female.GTEx)]
TCGA_TPM.male = TCGA_TPM[,which(colnames(TCGA_TPM) %in% male.TCGA)]
TCGA_TPM.female = TCGA_TPM[,which(colnames(TCGA_TPM) %in% female.TCGA)]

##### Set Y gene expression to "NA" in female samples
GTEx_TPM.female[which(rownames(GTEx_TPM.female) %in% Y_genes.GTEx),] = "NA"
TCGA_TPM.female[which(rownames(TCGA_TPM.female) %in% Y_genes.TCGA),] = "NA"

# Remove version numbers (after ".") from gene names
rownames(GTEx_TPM.male) = gsub("\\..*","",rownames(GTEx_TPM.male)) 
rownames(GTEx_TPM.female) = gsub("\\..*","",rownames(GTEx_TPM.female)) 
rownames(TCGA_TPM.male) = gsub("\\..*","",rownames(TCGA_TPM.male)) 
rownames(TCGA_TPM.female) = gsub("\\..*","",rownames(TCGA_TPM.female))

write.table(GTEx_TPM.male, "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/lioness_expression/GTEx_TPM_male.txt", row.names=T, col.names=T, quote=F)
write.table(GTEx_TPM.female, "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/lioness_expression/GTEx_TPM_female.txt", row.names=T, col.names=T, quote=F)
write.table(TCGA_TPM.male, "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/lioness_expression/TCGA_TPM_male.txt", row.names=T, col.names=T, quote=F)
write.table(TCGA_TPM.female, "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/lioness_expression/TCGA_TPM_female.txt", row.names=T, col.names=T, quote=F)

