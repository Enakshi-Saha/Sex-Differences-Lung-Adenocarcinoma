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

# Load GTEx and TCGA RSE object
GTEx_lung = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/GTEx_RSE.RData"))
TCGA_lung = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/TCGA_RSE.RData"))

# GTEx & TCGA raw counts
assays(GTEx_lung)$counts <- recount3::transform_counts(GTEx_lung)
assays(TCGA_lung)$counts <- recount3::transform_counts(TCGA_lung)

##### Check for batch effects in raw counts
GTEx_phenotypes = GTEx_lung@colData
# GTEx Batch variable: SMNABTCHT
GTEx_batch = GTEx_phenotypes$gtex.smnabtcht
# PCA to check batch effect
cols = factor(GTEx_batch)
library(ggplot2)
logCount = log2(1+assays(GTEx_lung)$counts)
pca.GTEx <- prcomp(t(logCount))
U <- pca.GTEx$x
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], main = "PCA of GTEx: log2(1+raw count)", colour = cols, xlab = "PC1", ylab = "PC2") + coord_equal()

##### Remove batch effects from raw counts
# GTEx Batch variable: SMNABTCHT
GTEx_batch = GTEx_lung@colData$gtex.smnabtcht
# Remove batch effect, controlling for sex
GTEx_raw = assays(GTEx_lung)$raw_counts
GTEx_gender = GTEx_lung@colData$gtex.sex
GTEx_raw <- ComBat_seq(GTEx_raw, batch=GTEx_batch, group=GTEx_gender, full_mod=TRUE)
logCount = log2(1+GTEx_raw)
pca.GTEx <- prcomp(t(logCount))
U <- pca.GTEx$x
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], main = "GTEx: batch-corrected log2(1+raw count)", colour = cols, xlab = "PC1", ylab = "PC2") + coord_equal()

#########################################
##### Check for batch effects in raw counts
TCGA_phenotypes = TCGA_lung@colData[which(TCGA_lung@colData$tcga.cgc_sample_sample_type == "Primary Tumor"),]
# TCGA Batch variable
TCGA_batch = TCGA_phenotypes$tcga.cgc_case_batch_number
# PCA to check batch effect
cols = factor(TCGA_batch)
library(ggplot2)
TCGA_raw_tumor = assays(TCGA_lung)$counts[,which(TCGA_lung@colData$tcga.cgc_sample_sample_type == "Primary Tumor")]
logCount = log2(1+TCGA_raw_tumor)
pca.TCGA <- prcomp(t(logCount))
U <- pca.TCGA$x
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], main = "PCA of TCGA: log2(1+raw count)", colour = cols, xlab = "PC1", ylab = "PC2") + coord_equal()

##### Remove batch effects from raw counts
TCGA_gender = TCGA_phenotypes$tcga.cgc_case_gender
TCGA_gender[is.na(TCGA_gender)] = "NA"
TCGA_raw <- ComBat_seq(TCGA_raw_tumor, batch=TCGA_batch, group=TCGA_gender, full_mod=TRUE)
logCount = log2(1+TCGA_raw)
pca.TCGA <- prcomp(t(logCount))
U <- pca.TCGA$x
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], main = "TCGA: batch-corrected log2(1+raw count)", colour = cols, xlab = "PC1", ylab = "PC2") + coord_equal()

