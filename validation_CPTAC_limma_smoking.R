#########################################
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
library(plyr)

# Load indegree
indegree = fread("/home/ubuntu/CPTAC2020_inDegree.txt")
head(indegree[,1:4])
genes = indegree$V1

# Get chromosome information
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
names(gene_info.TCGA)
head(gene_info.TCGA)
chr_loc = rep(0, length(genes))
for (i in 1:length(genes)){
  index = which(gene_info.TCGA$gene_name == genes[i])[1]
  chr_loc[i]=as.character(gene_info.TCGA$seqnames[index])
}
chr_loc[1:100]

# Check number of Y genes
which(chr_loc == "chrY")
indegree = indegree[,-1]
rownames(indegree) = genes

# Get clinical data
clinical = read.delim("~/validation_data/CPTAC2020/luad_cptac_2020/data_clinical_patient.txt", comment.char="#")
clinical$PATIENT_ID = unlist(lapply(strsplit(clinical$PATIENT_ID, split = "-"), function(x){paste(x, collapse = ".")}))
clinical = clinical[match(colnames(indegree),clinical$PATIENT_ID),]

gender = factor(clinical$SEX, levels = c("Male", "Female"))
age = clinical$AGE
age[which(is.na(age))] = median(age[-which(is.na(age))])
ethnicity = clinical$ETHNICITY
ethnicity[which(ethnicity %in% c("Caucasian", "European", "White"))] = "White"
ethnicity[which(ethnicity %in% c("Asian", "Han"))] = "Asian"
ethnicity[which(ethnicity %in% c("Black", "Hispanic"))] = "Others"
ethnicity[which(is.na(ethnicity))] = "Na"
ethnicity = factor(ethnicity, levels = c("White", "Asian", "Others", "Na"))
bmi = clinical$BMI
bmi[which(is.na(bmi))] = median(bmi[-which(is.na(bmi))])
tumor_stage = clinical$STAGE
tumor_stage[which(is.na(tumor_stage))] = "NA"
tumor_stage[which(tumor_stage %in% c("1","1A", "1B"))] = "stageI"
tumor_stage[which(tumor_stage %in% c("2A", "2B"))] = "stageII"
tumor_stage[which(tumor_stage %in% c("3","3A"))] = "stageIII"
tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII", "NA"))
smoking = clinical$SMOKING_STATUS
smoking[which(is.na(smoking))] = "NA"
smoking = factor(smoking, levels = c("Non-Smoker", "Smoker", "NA"))

#### Design Matrix
design = model.matrix(~ age + gender*smoking + ethnicity + bmi + tumor_stage + smoking)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)
rownames(fit$coefficients) = rownames(indegree)
#save(fit, file = "/home/esaha/limma_tables/TMM/unfiltered/netZooPy0.9.9/validation_CPTAC2020_limma.RData")
#save(design, file = "/home/esaha/limma_tables/TMM/unfiltered/netZooPy0.9.9/validation_CPTAC2020_design.RData")

#################

# Save table for sex difference: genderfemale
tb = topTable(fit,coef="genderFemale",number=Inf)

############# GSEA ###################
# Rank genes in limma table
indegree_rank <- setNames(object=tb[,"t"], rownames(tb))
head(indegree_rank)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
# Load reactome pathways
# pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
# Load GO pathways
# pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c5.go.bp.v2022.1.Hs.symbols.gmt")

fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
head(fgseaRes)

# Cancer Pathways
fgseaRes[which(fgseaRes$pathway %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_PATHWAYS_IN_CANCER", "KEGG_SMALL_CELL_LUNG_CANCER",
                                       "KEGG_NON_SMALL_CELL_LUNG_CANCER", "KEGG_WNT_SIGNALING_PATHWAY", "KEGG_ERBB_SIGNALING_PATHWAY",
                                       "KEGG_P53_SIGNALING_PATHWAY"))]


# Immune Pathways

rbind(fgseaRes[grep("immun", fgseaRes$pathway, ignore.case = T),], fgseaRes[grep("kine", fgseaRes$pathway, ignore.case = T),])

write.table(tb, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_CPTAC2020_nonsmoker.txt", col.names = T, row.names = T)
save(fgseaRes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_CPTAC2020_nonsmoker.RData")

