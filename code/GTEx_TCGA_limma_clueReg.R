# Must set seed
rm(list = ls())
set.seed(123)

# Load GTEx Male & Female Indegree

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("variancePartition")
#BiocManager::install("doParallel")

library(fgsea)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(recount)

##################################
# Load GTEx indegree and phenotypes
indegree_GTEx_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_male.txt"))
indegree_GTEx_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_female.txt"))
genes = indegree_GTEx_male$V1

indegree_GTEx = cbind(indegree_GTEx_male[,-1], indegree_GTEx_female[,-1])
rownames(indegree_GTEx) = genes

GTEx_IDs = colnames(indegree_GTEx)
# Format GTEx IDs to match with phenotypic data: GTEx
GTEx_IDs = unlist(lapply(strsplit(GTEx_IDs, split=".", fixed = T),
                         function(x){paste(x[1:2], collapse= "-")}))

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(indegree_GTEx), gene_info.TCGA$gene_id)]
chr_loc = gene_info.TCGA$seqnames[match(rownames(indegree_GTEx), gene_info.TCGA$gene_id)]


# Remove Y genes
indegree_GTEx = indegree_GTEx[-which(chr_loc == "chrY"),]
genes = genes[-which(chr_loc == "chrY")]
rownames(indegree_GTEx) = genes

# Get phenotypes: GTEx
GTEx_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_phenotypes.txt", sep=""))
GTEx_phenotypes = GTEx_phenotypes[match(GTEx_IDs, GTEx_phenotypes$SUBJID),]
GTEx_phenotypes = GTEx_phenotypes[,c("gtex.sex", "AGE", "RACE", "MHSMKSTS", "TRISCH", "gtex.smrin", "gtex.smnabtcht")]
colnames(GTEx_phenotypes) = c("gender", "age", "race", "smoking", "ischemic_time", "rna_degrad", "batch")

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("MALE", "FEMALE"))
gender_GTEx = gender

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)
age_GTEx = age

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"
race_GTEx = factor(race)

smoking_status = GTEx_phenotypes$smoking
smoking_status[-which(smoking_status %in% c("No", "Yes"))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))
smoking_GTEx = smoking_status

###############################
# Load TCGA indegree and phenotypes
# Load indegree matrices
indegree_TCGA_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_male.txt"))
indegree_TCGA_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_female.txt"))

genes = indegree_TCGA_male$V1

indegree_TCGA = cbind(indegree_TCGA_male[,-1], indegree_TCGA_female[,-1])
rownames(indegree_TCGA) = genes

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(indegree_TCGA), gene_info.TCGA$gene_id)]
chr_loc = gene_info.TCGA$seqnames[match(rownames(indegree_TCGA), gene_info.TCGA$gene_id)]

# Remove Y genes
indegree_TCGA = indegree_TCGA[-which(chr_loc == "chrY"),]
genes = genes[-which(chr_loc == "chrY")]
rownames(indegree_TCGA) = genes

# remove X from the starting of TCGA IDs
colnames(indegree_TCGA)[which(substring(colnames(indegree_TCGA), 1,1) == "X")] = substring(colnames(indegree_TCGA)[which(substring(colnames(indegree_TCGA), 1,1) == "X")], 2, length(colnames(indegree_TCGA)[which(substring(colnames(indegree_TCGA), 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(colnames(indegree_TCGA), split=".", fixed = T),
                         function(x){paste(x, collapse ="-")}))

# Get phenotypes: TCGA
TCGA_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_phenotypes.txt", sep=""))
TCGA_phenotypes = TCGA_phenotypes[match(TCGA_IDs, rownames(TCGA_phenotypes)),]
TCGA_phenotypes = TCGA_phenotypes[,c("tcga.gdc_cases.demographic.gender", "tcga.gdc_cases.demographic.race", 
                                     "tcga.xml_age_at_initial_pathologic_diagnosis", "tcga.gdc_cases.samples.sample_type", 
                                     "tcga.gdc_cases.diagnoses.tumor_stage", "tcga.xml_tobacco_smoking_history")]
colnames(TCGA_phenotypes) = c("gender", "race", "age", "sample_type", "tumor_stage", "smoking_status")
head(TCGA_phenotypes)
dim(TCGA_phenotypes)

# Define the covariates: gender, age, race, smoking, tumor_stage
gender = TCGA_phenotypes$gender
gender[which(gender == "male")] = "MALE" 
gender[which(gender == "female")] = "FEMALE" 
gender = factor(gender, levels = c("MALE", "FEMALE"))
gender_TCGA = gender

race = TCGA_phenotypes$race
race[which(race != "black or african american" & race != "white")] = "others"
race = factor(race)
race_TCGA = race

age <- as.numeric(as.character(TCGA_phenotypes$age))
age[which(is.na(age))] = mean(age,na.rm=TRUE)
age_TCGA = age

smoking_status = TCGA_phenotypes$smoking_status
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status[which(smoking_status == 1)] = "No"
smoking_status[which(smoking_status %in% 2:5)] = "Yes"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))
smoking_TCGA = smoking_status

###################################
###### Combined Indegree and phenotypes ######
indegree = cbind(indegree_GTEx,indegree_TCGA)
head(indegree[,1:4])

age = c(age_GTEx, age_TCGA)
gender = factor(c(gender_GTEx, gender_TCGA), levels = c("MALE", "FEMALE"))
race = c(race_GTEx, race_TCGA)
smoking_status = c(smoking_GTEx, smoking_TCGA)
disease = c(rep("healthy", ncol(indegree_GTEx)), rep("tumor", ncol(indegree_TCGA)))
disease = factor(disease, levels = c("healthy", "tumor"))

#### Design Matrix
design = model.matrix(~ age + gender + race + smoking_status +  gender*disease)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderFEMALE
tb = topTable(fit,coef="diseasetumor",number=Inf)
tb$chr = c()
tb$gene_name = c()
for (i in 1:nrow(tb)){
  index = match(rownames(tb)[i], gene_info.TCGA$gene_id)
  tb$chr[i]=gene_info.TCGA$seqnames[index]
  tb$gene_name[i]=gene_info.TCGA$gene_name[index]
}
head(tb)

pos_genes = tb$gene_name[which(tb$t > 0 & tb$adj.P.Val <0.05)]
neg_genes = tb$gene_name[which(tb$t < 0 & tb$adj.P.Val <0.05)]

write.table(pos_genes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/clueReg/positive_genes_male_0.05.txt", row.names = F, col.names = F, quote = F)
write.table(neg_genes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/clueReg/negative_genes_male_0.05.txt", row.names = F, col.names = F, quote = F)
