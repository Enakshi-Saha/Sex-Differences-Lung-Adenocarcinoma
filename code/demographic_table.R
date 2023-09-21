###### Column 1: GTEx ######
GTEx_phenotypes <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_phenotypes.txt", sep="")

nrow(GTEx_phenotypes)
table(GTEx_phenotypes$SEX)
table(GTEx_phenotypes$SEX)/nrow(GTEx_phenotypes)

mean(GTEx_phenotypes$AGE)
sd(GTEx_phenotypes$AGE)
range(GTEx_phenotypes$AGE)

table(GTEx_phenotypes$RACE)
table(GTEx_phenotypes$RACE)/nrow(GTEx_phenotypes)

table(GTEx_phenotypes$MHSMKSTS)
table(GTEx_phenotypes$MHSMKSTS)/nrow(GTEx_phenotypes)

###### Column 2: TCGA ######
TCGA_phenotypes <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_phenotypes.txt", sep="")

nrow(TCGA_phenotypes)
table(TCGA_phenotypes$tcga.cgc_case_gender)
table(TCGA_phenotypes$tcga.cgc_case_gender)/nrow(TCGA_phenotypes)

age = na.omit(as.numeric(as.character(TCGA_phenotypes$tcga.cgc_case_age_at_diagnosis)))
mean(age)
sd(age)
range(age)

race = TCGA_phenotypes$tcga.gdc_cases.demographic.race
table(race)
table(race)/nrow(TCGA_phenotypes)

smoking_status = TCGA_phenotypes$tcga.xml_tobacco_smoking_history
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status[which(smoking_status == 1)] = "No"
smoking_status[which(smoking_status %in% 2:5)] = "Yes"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))
table(smoking_status)
table(smoking_status)/nrow(TCGA_phenotypes)

tumor_stage = TCGA_phenotypes$tcga.gdc_cases.diagnoses.tumor_stage
tumor_stage[which(tumor_stage == "stage i" | tumor_stage == "stage ia" | tumor_stage == "stage ib")] = "stageI"
tumor_stage[which(tumor_stage == "stage ii" | tumor_stage == "stage iia" | tumor_stage == "stage iib")] = "stageII"
tumor_stage[which(tumor_stage == "Stage iii" | tumor_stage == "stage iiia" | tumor_stage == "stage iiib")] = "stageIII"
tumor_stage[which(tumor_stage == "stage iv")] = "stageIV"
tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII", "stageIV", "not reported"))
table(tumor_stage)
table(tumor_stage)/nrow(TCGA_phenotypes)

###### Column 4: GSE68465 ######
GSE68465_phenotypes <- read.csv("~/validation_data/validation_newdata2023/GSE68465_phenotypes.txt", sep="")

nrow(GSE68465_phenotypes)
table(GSE68465_phenotypes$Sex.ch1)
table(GSE68465_phenotypes$Sex.ch1)/nrow(GSE68465_phenotypes)

mean(GSE68465_phenotypes$age.ch1)
sd(GSE68465_phenotypes$age.ch1)
range(GSE68465_phenotypes$age.ch1)

table(GSE68465_phenotypes$race.ch1)
table(GSE68465_phenotypes$race.ch1)/nrow(GSE68465_phenotypes)

smoking = GSE68465_phenotypes$smoking_history.ch1
smoking[which(smoking %in% c("Currently smoking", "Smoked in the past"))] = "Ever"
smoking[which(smoking == "Never smoked")] = "Never"
smoking[-which(smoking %in% c("Ever", "Never"))] = "Unknown"
table(smoking)
table(smoking)/nrow(GSE68465_phenotypes)

tumor_stage = GSE68465_phenotypes$disease_stage.ch1
tumor_stage = substr(stage, nchar(stage), nchar(stage))
table(tumor_stage)
table(tumor_stage)/nrow(GSE68465_phenotypes)
