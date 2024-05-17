# Must set seed
rm(list = ls())
set.seed(123)
library(ggplot2)
library(data.table)
library(recount)
library(gghalves)

# TCGA: XIST distribution by sex and smoking status
# Load indegree matrices
indegree_TCGA_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_male.txt"))
indegree_TCGA_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_female.txt"))

genes = indegree_TCGA_male$V1

indegree_TCGA = cbind(indegree_TCGA_male[,-1], indegree_TCGA_female[,-1])
rownames(indegree_TCGA) = genes

# Load expression matrices
TCGA_TPM = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_TPM.txt"))
TCGA_TPM$V1 = gsub("\\..*","", TCGA_TPM$V1)
TCGA_TPM = TCGA_TPM[match(genes, TCGA_TPM$V1),match(colnames(indegree_TCGA), colnames(TCGA_TPM))]
rownames(TCGA_TPM) = genes

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(indegree_TCGA), gene_info.TCGA$gene_id)]

xist_indegree = as.numeric(indegree_TCGA[which(gene_name == "XIST"),])
xist_expression = as.numeric(TCGA_TPM[which(gene_name == "XIST"),])

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

race = TCGA_phenotypes$race
race[which(race != "black or african american" & race != "white")] = "others"
race = factor(race)

age <- as.numeric(as.character(TCGA_phenotypes$age))
age[which(is.na(age))] = mean(age,na.rm=TRUE)

smoking_status = TCGA_phenotypes$smoking_status
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status[which(smoking_status == 1)] = "No"
smoking_status[which(smoking_status %in% 2:5)] = "Yes"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

tumor_stage = TCGA_phenotypes$tumor_stage
tumor_stage[which(tumor_stage == "stage i" | tumor_stage == "stage ia" | tumor_stage == "stage ib")] = "stageI"
tumor_stage[which(tumor_stage == "stage ii" | tumor_stage == "stage iia" | tumor_stage == "stage iib")] = "stageII"
tumor_stage[which(tumor_stage == "Stage iii" | tumor_stage == "stage iiia" | tumor_stage == "stage iiib")] = "stageIII"
tumor_stage[which(tumor_stage == "stage iv")] = "stageIV"
#tumor_stage[which(tumor_stage %in% c("stageIII", "stageIV"))] = "stageIII-IV"
#tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII-IV", "not reported"))
tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII", "stageIV", "not reported"))

# Create dataframe for boxplot
fulldata = data.frame(xist_indegree)
fulldata$expr = xist_expression
fulldata$gender = factor(gender, levels = c("FEMALE", "MALE"))
fulldata$smoking_status = smoking_status

# Violinplot with jitters
maintitle = "TCGA: Indegree Distribution of XIST"
ggplot(fulldata, aes(x = smoking_status, y = xist_indegree, fill = gender)) + geom_violin() + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("indegree")
maintitle = "TCGA: Expression Distribution of XIST"
ggplot(fulldata, aes(x = smoking_status, y = expr, col = gender)) + geom_violin() + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("expression")

# Boxplot
maintitle = "TCGA: Indegree Distribution of XIST"
ggplot(fulldata, aes(x = smoking_status, y = xist_indegree, fill = gender)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("indegree")
maintitle = "TCGA: Expression Distribution of XIST"
ggplot(fulldata, aes(x = smoking_status, y = expr, col = gender)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("expression")

# Wilcoxon test
sample1 = fulldata$expr[which(fulldata$gender == "MALE" & fulldata$smoking_status == "No")]
sample2 = fulldata$expr[which(fulldata$gender == "FEMALE" & fulldata$smoking_status == "No")]
wilcox.test(sample1, sample2)

sample1 = fulldata$expr[which(fulldata$gender == "MALE" & fulldata$smoking_status == "Yes")]
sample2 = fulldata$expr[which(fulldata$gender == "FEMALE" & fulldata$smoking_status == "Yes")]
wilcox.test(sample1, sample2)

