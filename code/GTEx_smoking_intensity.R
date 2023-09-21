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
GTEx_phenotypes = GTEx_phenotypes[,c("gtex.sex", "AGE", "RACE", "MHSMKSTS", "TRISCH", "MHSMKYRS", "MHSMKNMB", "MHSMKPRD")]
colnames(GTEx_phenotypes) = c("gender", "age", "race", "smoking", "ischemic_time", "smoke_years", "smoke_number", "smoke_period")

gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
GTEx_gender = factor(gender, levels = c("MALE", "FEMALE"))

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)
GTEx_age = age

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"
GTEx_race <- as.factor(race)

smoking_status = GTEx_phenotypes$smoking
smoking_status[which(is.na(smoking_status))] = "Unknown"
GTEx_smoking_status = smoking_status

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})

smoke_years = GTEx_phenotypes$smoke_years
smoke_num = GTEx_phenotypes$smoke_number
smoke_period = factor(GTEx_phenotypes$smoke_period)
total_smoke = rep("NA", length(gender))
total_smoke[smoke_period == "Day"] = 365*smoke_years[smoke_period == "Day"]*smoke_num[smoke_period == "Day"]
total_smoke[smoke_period == "Week"] = 52*smoke_years[smoke_period == "Week"]*smoke_num[smoke_period == "Week"]
total_smoke[smoke_period == "Month"] = 12*smoke_years[smoke_period == "Month"]*smoke_num[smoke_period == "Month"]
total_smoke = log10(as.numeric(as.character(total_smoke)))

##### Get pathway score
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
indegree = indegree_GTEx
pathway_score = c()
for (i in 1:length(pathways)){
  indegrees = indegree[which(gene_name %in% pathways[i][[1]]),]
  pathway_score = rbind(pathway_score,unlist(apply(indegrees, MARGIN = 2, FUN = mean)))
}
rownames(pathway_score) = names(pathways)

cor.test(smoke_years, pathway_score[which(rownames(pathway_score) == "KEGG_NON_SMALL_CELL_LUNG_CANCER"),], use = "pairwise.complete.obs")
cor.test(smoke_years, pathway_score[which(rownames(pathway_score) == "KEGG_WNT_SIGNALING_PATHWAY"),], use = "pairwise.complete.obs")

library(ggplot2)
fulldata = data.frame(as.numeric(total_smoke))
fulldata$gender = factor(gender)
ggplot(fulldata, aes(x=gender, y=total_smoke, fill = gender)) + geom_boxplot() + ggtitle("GTEx: Smoking History by Sex") + xlab("") + ylab("log10(lifetime number of cigarettes)")

fulldata$score = pathway_score[which(rownames(pathway_score) == "KEGG_NON_SMALL_CELL_LUNG_CANCER"),]
rownames(fulldata) = colnames(pathway_score)
fulldata = na.omit(fulldata)
fulldata = fulldata[-which(fulldata$as.numeric.total_smoke. == -Inf),]
wilcox.test(fulldata$as.numeric.total_smoke.[fulldata$gender == "Male"], fulldata$as.numeric.total_smoke.[fulldata$gender == "Female"], alternative = "greater")

cor.test(fulldata$as.numeric.total_smoke., fulldata$score)
cor.test(fulldata$as.numeric.total_smoke[which(fulldata$gender == "MALE")], fulldata$score[which(fulldata$gender == "MALE")])
cor.test(fulldata$as.numeric.total_smoke[which(fulldata$gender == "FEMALE")], fulldata$score[which(fulldata$gender == "FEMALE")])

############################
# Smoking intensity adjusted GSEA
fulldata = data.frame(as.numeric(total_smoke))
rownames(fulldata) = colnames(pathway_score)
fulldata = na.omit(fulldata)
included_samples = rownames(fulldata)[-which(fulldata$as.numeric.total_smoke. == -Inf)]


indegree0 = indegree[,which(colnames(indegree) %in% included_samples)]
gender0 = gender[which(colnames(indegree) %in% included_samples)]
gender0 = factor(gender0, levels = c("MALE", "FEMALE"))
age0 = age[which(colnames(indegree) %in% included_samples)]
race0 = race[which(colnames(indegree) %in% included_samples)]
ischemic_time0 = ischemic_timeH[which(colnames(indegree) %in% included_samples)]
total_smoke0 = as.numeric(total_smoke[which(colnames(indegree) %in% included_samples)])

#### Design Matrix
design = model.matrix(~ age0 + gender0 + race0 + total_smoke0 + ischemic_time0 )
fit <- lmFit(indegree0, design)
fit <- eBayes(fit)

# Save table for sex difference: genderFEMALE
tb = topTable(fit,coef="gender0FEMALE",number=Inf)
tb$chr = c()
tb$gene_name = c()
for (i in 1:nrow(tb)){
  index = match(rownames(tb)[i], gene_info.TCGA$gene_id)
  tb$chr[i]=gene_info.TCGA$seqnames[index]
  tb$gene_name[i]=gene_info.TCGA$gene_name[index]
}
head(tb)

############# GSEA ###################
# Rank genes in limma table
indegree_rank <- setNames(object=tb[,"t"], tb$gene_name)
head(indegree_rank)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
head(fgseaRes)

# Cancer Pathways
fgseaRes[which(fgseaRes$pathway %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_PATHWAYS_IN_CANCER", "KEGG_SMALL_CELL_LUNG_CANCER",
                                       "KEGG_NON_SMALL_CELL_LUNG_CANCER", "KEGG_WNT_SIGNALING_PATHWAY", "KEGG_ERBB_SIGNALING_PATHWAY",
                                       "KEGG_P53_SIGNALING_PATHWAY"))]


