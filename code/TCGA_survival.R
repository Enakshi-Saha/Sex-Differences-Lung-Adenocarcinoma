# Survival Analysis with Pathway Score
set.seed(0)
rm(list = ls())
library(fgsea)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(recount3)
library(recount)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

# Load indegree matrices
indegree_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_male.txt"))
indegree_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_female.txt"))

genes = indegree_male$V1
indegree = cbind(indegree_male[,-1], indegree_female[,-1])
rownames(indegree) = genes

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(indegree), gene_info.TCGA$gene_id)]
chr_loc = gene_info.TCGA$seqnames[match(rownames(indegree), gene_info.TCGA$gene_id)]

# Remove Y genes
gene_name = gene_name[-which(chr_loc == "chrY")]
indegree = indegree[-which(chr_loc == "chrY"),]
genes = genes[-which(chr_loc == "chrY")]
rownames(indegree) = genes

##### Get pathway score
pathway_score = c()
for (i in 1:length(pathways)){
  indegrees = indegree[which(gene_name %in% pathways[i][[1]]),]
  pathway_score = rbind(pathway_score,unlist(apply(indegrees, MARGIN = 2, FUN = mean)))
}
rownames(pathway_score) = names(pathways)

# remove X from the starting of TCGA IDs
colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")] = substring(colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")], 2, length(colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(colnames(indegree), split=".", fixed = T),
                         function(x){paste(x, collapse ="-")}))

# Get phenotypes: TCGA
TCGA_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_phenotypes.txt", sep=""))
TCGA_phenotypes = TCGA_phenotypes[match(TCGA_IDs, rownames(TCGA_phenotypes)),]
TCGA_phenotypes = TCGA_phenotypes[,c("tcga.gdc_cases.demographic.gender", "tcga.gdc_cases.demographic.race", 
                                     "tcga.xml_age_at_initial_pathologic_diagnosis", "tcga.gdc_cases.samples.sample_type", 
                                     "tcga.gdc_cases.diagnoses.tumor_stage", "tcga.xml_tobacco_smoking_history",
                                     "tcga.xml_days_to_death", "tcga.xml_days_to_last_followup",
                                     "tcga.xml_vital_status", "tcga.cgc_drug_therapy_pharmaceutical_therapy_type")]
colnames(TCGA_phenotypes) = c("gender", "race", "age", "sample_type", "tumor_stage", "smoking_status", "death_days", "last_contact_days", "vital_status", "therapy")
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

therapy = TCGA_phenotypes$therapy
therapy[which(is.na(therapy))] = "NA"
therapy[-which(therapy %in% c("Chemotherapy", "NA"))] = "other"

death_days <- as.numeric(as.character(TCGA_phenotypes$death_days))
#death_days[which(is.na(death_days))] <- mean(death_days,na.rm=TRUE)

last_contact_days <- as.numeric(as.character(TCGA_phenotypes$last_contact_days))
#last_contact_days[which(is.na(last_contact_days))] <- mean(last_contact_days,na.rm=TRUE)

vital_status <- TCGA_phenotypes$vital_status
vital_status = factor(vital_status)

# Survival Analysis: Cox Regression

time <-  ifelse(is.na(death_days), yes=last_contact_days, no=death_days)
time <- time/365
status <- vital_status
status <- sub("Alive", 0, status)
status <- sub("Dead", 1, status)
status <- sub("[Discrepancy]", NA, status)
status <- as.numeric(status)

## Load survival package
library(survival)
SurvObj = Surv(time, status)
diff = survdiff(Surv(time,status)~gender)
pval <- round(1 - pchisq(diff$chisq, length(diff$n) - 1),3)
plot(survfit(Surv(time,status)~gender), main = "Plot of Survival Curves by Sex", xlab = "Length of Survival",ylab="Proportion of Individuals who have Survived",col=c("blue","red"))
legend("topright", legend=levels(gender),fill=c("blue","red"),bty="n")

# Fit Cox regression for every pathway
pvals_score = rep(0, nrow(pathway_score))
pvals_interaction = pvals_score
for (i in 1:nrow(pathway_score)){
  score = pathway_score[i,]
  cox_tcga <- coxph(SurvObj ~ age + gender + race + tumor_stage + smoking_status + therapy 
                    + gender*score)
  coeffs = summary(cox_tcga)$coefficients
  pvals_score[i] = coeffs["score",5]
  pvals_interaction[i] = coeffs["genderFEMALE:score",5]
}
names(pvals_score) = rownames(pathway_score)
names(pvals_interaction) = rownames(pathway_score)

which(pvals_interaction < 0.05)
which(pvals_score < 0.05)

########## Chemotherapy samples only ###########
# Fit Cox regression for every pathway
pvals_score = rep(0, nrow(pathway_score))
pvals_interaction = pvals_score
for (i in 1:nrow(pathway_score)){
  score = pathway_score[i,]
  cox_tcga <- coxph(SurvObj[which(therapy == "Chemotherapy")] ~ age[which(therapy == "Chemotherapy")] + gender[which(therapy == "Chemotherapy")] + race[which(therapy == "Chemotherapy")] + tumor_stage[which(therapy == "Chemotherapy")] + smoking_status[which(therapy == "Chemotherapy")] + 
              gender[which(therapy == "Chemotherapy")]*score[which(therapy == "Chemotherapy")])
  coeffs = summary(cox_tcga)$coefficients
  pvals_score[i] = coeffs[11,5]
  pvals_interaction[i] = coeffs[12,5]
}
names(pvals_score) = rownames(pathway_score)
names(pvals_interaction) = rownames(pathway_score)

which(pvals_interaction < 0.05)
which(pvals_score < 0.05)

