library(table1)

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform fisher test
    p <- fisher.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

########################################
# Table for GTEx
library(fgsea)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(recount)

# Load indegree matrices
indegree_GTEx_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_male.txt"))
indegree_GTEx_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_female.txt"))
genes = indegree_GTEx_male$V1

indegree_GTEx = cbind(indegree_GTEx_male[,-1], indegree_GTEx_female[,-1])
rownames(indegree_GTEx) = genes

GTEx_IDs = colnames(indegree_GTEx)
# Format GTEx IDs to match with phenotypic data: GTEx
GTEx_IDs = unlist(lapply(strsplit(GTEx_IDs, split=".", fixed = T),
                         function(x){paste(x[1:2], collapse= "-")}))

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
sex = gender

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

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})

dataset = data.frame(age, sex, race, smoking_status, ischemic_timeH )
head(dataset)

table1(~ age + race + smoking_status + ischemic_timeH | sex, data=dataset,  overall=F, extra.col=list(`P-value`=pvalue), topclass="Rtable1-grid Rtable1-shade Rtable1-times")

#####################################
# Table for TCGA
# Load indegree matrices
indegree_TCGA_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_male.txt"))
indegree_TCGA_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_female.txt"))

genes = indegree_TCGA_male$V1

indegree_TCGA = cbind(indegree_TCGA_male[,-1], indegree_TCGA_female[,-1])
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
sex = gender

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
tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII", "stageIV", "not reported"))

dataset = data.frame(age, sex, race, smoking_status, tumor_stage)
head(dataset)

table1(~ age + race + smoking_status + tumor_stage | sex, data=dataset,  overall=F, extra.col=list(`P-value`=pvalue), topclass="Rtable1-grid Rtable1-shade Rtable1-times")


#####################################
# Table for LGRC
# Load indegree
indegree = fread("/home/esaha/Lung_lioness/validation_networks/GSE47460/GSE47460_inDegree.txt")
head(indegree[,1:4])
genes = indegree$V1

# Load phenotypes
phenotypes = read.csv("/home/esaha/validation_data/GSE47460/GSE47460_phenotypes.txt", sep="")
phenotypes = phenotypes[match(colnames(indegree)[-1], rownames(phenotypes)),]

age = phenotypes$age
gender = factor(phenotypes$sex, levels = c("Male", "Female"))
sex = gender
smoking = phenotypes$smoking_status
smoking[is.na(smoking)] = "NA"
smoking = factor(smoking, levels = c("Never", "Ever", "NA"))

dataset = data.frame(age, sex, smoking)
head(dataset)

table1(~ age + smoking | sex, data=dataset,  overall=F, extra.col=list(`P-value`=pvalue), topclass="Rtable1-grid Rtable1-shade Rtable1-times")

#####################################
# Table for GSE68465
# Load indegree
indegree_male = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_male.txt")
indegree_female = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_female.txt")
indegree = cbind(indegree_male, indegree_female[,-1])
head(indegree[,1:4])
genes = indegree$V1

# Get clinical data
clinical = read.csv("~/validation_data/validation_newdata2023/GSE68465_phenotypes.txt", sep="")
clinical = clinical[match(colnames(indegree)[-1],clinical$geo_accession),]
head(clinical[,1:5])

gender = clinical$Sex.ch1
gender = factor(gender, levels = c("Male", "Female"))
sex = gender

race = clinical$race.ch1
race[which(race == "Not Reported")] = "Unknown"
race[-which(race %in% c("White", "Black or African American", "Unknown"))] = "others"
race = factor(race, levels = c("White", "Black or African American", "others", "Unknown"))

age = clinical$age.ch1

smoking = clinical$smoking_history.ch1
smoking[which(smoking %in% c("Currently smoking", "Smoked in the past"))] = "Ever"
smoking[which(smoking == "Never smoked")] = "Never"
smoking[-which(smoking %in% c("Ever", "Never"))] = "Unknown"
smoking = factor(smoking, levels = c("Never", "Ever", "Unknown"))
table(smoking)

stage = clinical$disease_stage.ch1
tumor_stage = substr(stage, nchar(stage), nchar(stage))
tumor_stage[which(tumor_stage == "p")] = "NA"
tumor_stage = factor(tumor_stage, levels = c("1", "2", "3", "4", "NA"))

#### Design Matrix
dataset = data.frame(age, sex, tumor_stage, race, smoking)
head(dataset)

table1(~ age + tumor_stage + race + smoking | sex, data=dataset,  overall=F, extra.col=list(`P-value`=pvalue), topclass="Rtable1-grid Rtable1-shade Rtable1-times")


