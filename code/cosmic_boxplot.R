set.seed(0)
rm(list = ls())

library(fgsea)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)

# Get oncogenes & tumor suppressor genes
cancer_genes <- read.delim("/home/ubuntu/cancer_genes_cosmic.txt")
head(cancer_genes)
oncogenes = cancer_genes$geneNames[which(cancer_genes$gene_type == "oncogene")]
tsg = cancer_genes$geneNames[which(cancer_genes$gene_type == "TSG")]

# Get limma matrices
GTEx_nonsmoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_GTEx_nonsmoker.txt", sep="")
GTEx_smoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_GTEx_smoker.txt", sep="")
TCGA_nonsmoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_TCGA_nonsmoker.txt", sep="")
TCGA_smoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_TCGA_smoker.txt", sep="")

GTEx_smoker = GTEx_smoker[match(GTEx_nonsmoker$gene_name,GTEx_smoker$gene_name),]
TCGA_nonsmoker = TCGA_nonsmoker[match(GTEx_nonsmoker$gene_name,TCGA_nonsmoker$gene_name),]
TCGA_smoker = TCGA_smoker[match(GTEx_nonsmoker$gene_name,TCGA_smoker$gene_name),]

tab_nonsmoker = data.frame(cbind(GTEx_nonsmoker$t, TCGA_nonsmoker$t))
colnames(tab_nonsmoker) = c("GTEx_nonsmoker", "TCGA_nonsmoker")
tab_nonsmoker$gene_name = GTEx_nonsmoker$gene_name
head(tab_nonsmoker)

tab_smoker = data.frame(cbind(GTEx_smoker$t, TCGA_smoker$t))
colnames(tab_smoker) = c("GTEx_smoker", "TCGA_smoker")
tab_smoker$gene_name = GTEx_nonsmoker$gene_name
head(tab_smoker)

gene_type = "Oncogene"
category = "Nonsmokers"
tab = tab_nonsmoker
tab_subset = tab[which(tab$gene_name %in% oncogenes),]
fulldata = data.frame(c(tab_subset[,1], tab_subset[,2]))
colnames(fulldata) = "tStat"
fulldata$sample_type = c(rep("Healthy", nrow(tab_subset)), rep("Tumor", nrow(tab_subset)))

maintitle = paste(gene_type, "Targeting:", category, "(Female - Male)", sep = " ")
ggplot(fulldata, aes(x = sample_type, y = tStat)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("gene t-statistic (female - male)") + geom_hline(yintercept=0, linetype="dashed",color = "red")
# Violinplot with jitters
ggplot(fulldata, aes(x = sample_type, y = tStat)) + geom_violin() + geom_jitter() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("gene t-statistic (female - male)") + geom_hline(yintercept=0, linetype="dashed",color = "red")


gene_type = "Tumor Suppresor Gene"
category = "Smokers"
tab = tab_smoker
tab_subset = tab[which(tab$gene_name %in% tsg),]
fulldata = data.frame(c(tab_subset[,1], tab_subset[,2]))
colnames(fulldata) = "tStat"
fulldata$sample_type = c(rep("Healthy", nrow(tab_subset)), rep("Tumor", nrow(tab_subset)))

maintitle = paste(gene_type, "Targeting:", category, "(Female - Male)", sep = " ")
ggplot(fulldata, aes(x = sample_type, y = tStat)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("gene t-statistic (female - male)") + geom_hline(yintercept=0, linetype="dashed",color = "red")
# Violinplot with jitters
ggplot(fulldata, aes(x = sample_type, y = tStat)) + geom_violin() + geom_jitter() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("gene t-statistic (female - male)") + geom_hline(yintercept=0, linetype="dashed",color = "red")
