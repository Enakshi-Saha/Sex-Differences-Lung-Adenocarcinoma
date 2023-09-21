set.seed(0)
rm(list = ls())

library(fgsea)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)

drug_genes = c("EGFR", "ALK", "MET", "ERBB2", "ERBB3", "ROS", "RET", 
               "FGFR1", "DDR2", "IGF1R", "KRAS", "NF1", "HRAS", "NRAS",
               "RASA1", "BRAF", "PIK3CA", "PTEN", "AKT1", "AKT2", "AKT3",
               "TSC1", "TSC2", "LKB1", "TP53", "MDM2", "CDKNA2", "RB1", "CCNE1",
               "MYC", "MYCN", "MYCL", "KRAS")

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

# Isolate drug genes
indegree_TCGA = indegree_TCGA[which(gene_name %in% drug_genes),]
genes = genes[which(gene_name %in% drug_genes)]
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

race = TCGA_phenotypes$race
race[which(race != "black or african american" & race != "white")] = "others"
race = factor(race)

age <- as.numeric(as.character(TCGA_phenotypes$age))
age[which(is.na(age))] = mean(age,na.rm=TRUE)

smoking_status = TCGA_phenotypes$smoking_status
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status[which(smoking_status == 1)] = "No"
smoking_status[which(smoking_status %in% 2:5)] = "Yes"
#smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))
smoking_status = factor(smoking_status, levels = c("Yes", "No", "Unknown"))

tumor_stage = TCGA_phenotypes$tumor_stage
tumor_stage[which(tumor_stage == "stage i" | tumor_stage == "stage ia" | tumor_stage == "stage ib")] = "stageI"
tumor_stage[which(tumor_stage == "stage ii" | tumor_stage == "stage iia" | tumor_stage == "stage iib")] = "stageII"
tumor_stage[which(tumor_stage == "Stage iii" | tumor_stage == "stage iiia" | tumor_stage == "stage iiib")] = "stageIII"
tumor_stage[which(tumor_stage == "stage iv")] = "stageIV"
#tumor_stage[which(tumor_stage %in% c("stageIII", "stageIV"))] = "stageIII-IV"
#tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII-IV", "not reported"))
tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII", "stageIV", "not reported"))

design = model.matrix(~ age + gender + race + gender*smoking_status + tumor_stage)
indegree = indegree_TCGA
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderFEMALE
tb = topTable(fit,coef="genderFEMALE",number=Inf)
tb$chr = c()
tb$gene_name = c()
for (i in 1:nrow(tb)){
  index = match(rownames(tb)[i], gene_info.TCGA$gene_id)
  tb$chr[i]=gene_info.TCGA$seqnames[index]
  tb$gene_name[i]=gene_info.TCGA$gene_name[index]
}
head(tb)
tb_TCGA = tb

################################
# Load indegree
indegree_male = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_male.txt")
indegree_female = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_female.txt")
indegree = cbind(indegree_male, indegree_female[,-1])
head(indegree[,1:4])
genes = indegree$V1

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(genes, gene_info.TCGA$gene_name)]
chr_loc = gene_info.TCGA$seqnames[match(genes, gene_info.TCGA$gene_name)]

# Remove Y genes
indegree = indegree[which(gene_name %in% drug_genes),]
genes = genes[which(gene_name %in% drug_genes)]
indegree = data.frame(indegree[,-1])
rownames(indegree) = genes

# Get clinical data
clinical = read.csv("~/validation_data/validation_newdata2023/GSE68465_phenotypes.txt", sep="")
clinical = clinical[match(colnames(indegree),clinical$geo_accession),]
head(clinical[,1:5])

gender = clinical$Sex.ch1
gender = factor(gender, levels = c("Male", "Female"))

race = clinical$race.ch1
race[which(race == "Not Reported")] = "Unknown"
race[-which(race %in% c("White", "Black or African American", "Unknown"))] = "others"
race = factor(race)

age = clinical$age.ch1

smoking = clinical$smoking_history.ch1
smoking[which(smoking %in% c("Currently smoking", "Smoked in the past"))] = "Ever"
smoking[which(smoking == "Never smoked")] = "Never"
smoking[-which(smoking %in% c("Ever", "Never"))] = "Unknown"
#smoking = factor(smoking, levels = c("Never", "Ever", "Unknown"))
smoking = factor(smoking, levels = c("Ever", "Never", "Unknown"))

stage = clinical$disease_stage.ch1
tumor_stage = substr(stage, nchar(stage), nchar(stage))
tumor_stage[which(tumor_stage %in% c("3", "4"))] = "3-4"
tumor_stage = factor(tumor_stage, levels = c("1", "2", "3-4", "p"))

#### Design Matrix
design = model.matrix(~ age + gender + tumor_stage + race + gender*smoking)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderfemale
tb_GSE68465 = topTable(fit,coef="genderFemale",number=Inf)

##########################################
tb_GSE68465 = tb_GSE68465[match(rownames(tb_GSE68465), tb_TCGA$gene_name),]

tb_subset = data.frame(cbind(tb_TCGA$t, tb_GSE68465$t, tb_TCGA$P.Value, tb_GSE68465$P.Value))
rownames(tb_subset) = tb_TCGA$gene_name
colnames(tb_subset) = c("TCGA", "GSE68465", "TCGA_pval", "GSE68465_pval")
tb_subset

#write.table(tb_subset, file = "/home/esaha/paper_plots/plots2023/drugGene_heatmap/TCGA_GSE68465_nonsmoker.txt", row.names = T, col.names = T)
write.table(tb_subset, file = "/home/esaha/paper_plots/plots2023/drugGene_heatmap/TCGA_GSE68465_smoker.txt", row.names = T, col.names = T)

###########################################

nonsmoker <- read.csv("~/paper_plots/plots2023/drugGene_heatmap/TCGA_GSE68465_nonsmoker.txt", sep="")
smoker <- read.csv("~/paper_plots/plots2023/drugGene_heatmap/TCGA_GSE68465_smoker.txt", sep="")

smoker = smoker[match(rownames(nonsmoker), rownames(smoker)),]

sig_genes = rownames(nonsmoker)[which(nonsmoker$TCGA_pval<0.05 | smoker$TCGA_pval<0.05)]

tab = cbind(nonsmoker$TCGA, nonsmoker$GSE68465, smoker$TCGA, smoker$GSE68465)[match(sig_genes, rownames(nonsmoker)),]
colnames(tab) = c("TCGA nonsmoker", "GSE68465 nonsmoker", "TCGA smoker", "GSE68465 smoker")
rownames(tab) = sig_genes

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/drugGene_heatmap/TCGA_GSE68465_drugGene_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 90, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

