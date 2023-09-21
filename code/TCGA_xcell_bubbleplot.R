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

# Load immune cell composition
TCGA_xcell <- read.csv("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/xcell/TCGA_xcell.txt", sep="")
rownames(TCGA_xcell) = TCGA_xcell$cell_type
TCGA_xcell = TCGA_xcell[,-1]

sample.TCGA = colnames(TCGA_xcell)
# remove X from the starting of TCGA IDs
sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")] = substring(sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")], 2, length(sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(sample.TCGA, split=".", fixed = T),
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

# Linear model coefficient of each cell type
TCGA_lm = apply(TCGA_xcell, MARGIN = 1, function(x){summary(lm(x ~ age + gender + race + tumor_stage + smoking_status))$coefficients[3,c(1,4)]})
rownames(TCGA_lm) = c("sex_coefficient", "p-value")
TCGA_lm = t(TCGA_lm)
TCGA_lm

# Bubbleplot plot by group
# importing the ggplot2 library
library(ggplot2)
# creating the dataframe from the above columns
fulldata = data.frame(TCGA_lm)
fulldata = fulldata[which(fulldata$p.value < 0.05),]
fulldata$cell_type = factor(rownames(fulldata))
fulldata$trend = rep("male", nrow(fulldata))
fulldata$trend[which(fulldata$sex_coefficient >= 0)] = "female"
fulldata$trend = factor(fulldata$trend, levels = c("female", "male"))
head(fulldata)

# Plot bubble plot
ggplot(fulldata, aes(x=p.value, y=cell_type, size=sex_coefficient, color = trend))+
  geom_point(alpha=0.8) + labs(title="TCGA: Sex Difference in Cell Composition",y="cell type", x = "pvalue")

