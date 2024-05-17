# Must set seed
rm(list = ls())
set.seed(123)
library(ggplot2)
library(data.table)
library(recount)

# GTEx: XIST distribution by sex and smoking status
# Load indegree matrices
indegree_GTEx_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_male.txt"))
indegree_GTEx_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_female.txt"))
genes = indegree_GTEx_male$V1

indegree_GTEx = cbind(indegree_GTEx_male[,-1], indegree_GTEx_female[,-1])
rownames(indegree_GTEx) = genes

# Load expression matrices
GTEx_TPM = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_TPM.txt"))
GTEx_TPM$V1 = gsub("\\..*","", GTEx_TPM$V1)
GTEx_TPM = GTEx_TPM[match(genes, GTEx_TPM$V1),]
GTEx_TPM = GTEx_TPM[,-1]
colnames(GTEx_TPM) = unlist(lapply(colnames(GTEx_TPM),function(x){substr(x,1,nchar(x)-2)}))
GTEx_TPM = GTEx_TPM[,match(colnames(indegree_GTEx), colnames(GTEx_TPM))]
rownames(GTEx_TPM) = genes

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

xist_indegree = as.numeric(indegree_GTEx[which(gene_name == "XIST"),])
xist_expression = as.numeric(GTEx_TPM[which(gene_name == "XIST"),])

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

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"

smoking_status = GTEx_phenotypes$smoking
smoking_status[-which(smoking_status %in% c("No", "Yes"))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

batch = factor(GTEx_phenotypes$batch)

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})

rin = as.numeric(as.character(GTEx_phenotypes$rna_degrad))

# Create dataframe for boxplot
fulldata = data.frame(xist_indegree)
fulldata$expr = xist_expression
fulldata$gender = factor(gender, levels = c("FEMALE", "MALE"))
fulldata$smoking_status = smoking_status

# Violinplot with jitters
maintitle = "GTEx: Indegree Distribution of XIST"
ggplot(fulldata, aes(x = smoking_status, y = xist_indegree, fill = gender)) + geom_violin() + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("indegree")
maintitle = "GTEx: Expression Distribution of XIST"
ggplot(fulldata, aes(x = smoking_status, y = expr, col = gender)) + geom_violin() + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("expression")

# Boxplot
ggplot(fulldata, aes(x = smoking_status, y = expr, fill = gender)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("indegree")

# Wilcoxon test
sample1 = fulldata$expr[which(fulldata$gender == "MALE" & fulldata$smoking_status == "No")]
sample2 = fulldata$expr[which(fulldata$gender == "FEMALE" & fulldata$smoking_status == "No")]
wilcox.test(sample1, sample2)

sample1 = fulldata$expr[which(fulldata$gender == "MALE" & fulldata$smoking_status == "Yes")]
sample2 = fulldata$expr[which(fulldata$gender == "FEMALE" & fulldata$smoking_status == "Yes")]
wilcox.test(sample1, sample2)

