# Examine leading genes in the Ribosome pathway
# Ribosome is male-biased in GTEx nonsmoker and female-biased in GTEx smoker and all TCGA
# Must set seed
rm(list = ls())
set.seed(123)
library(ggplot2)
library(data.table)
library(recount)
library(fgsea)
library(gghalves)

pathway = "KEGG_RIBOSOME"

GTEx_nonsmoker = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_nonsmoker.RData"))
GTEx_smoker = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_smoker.RData"))

TCGA_nonsmoker = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_nonsmoker.RData"))
TCGA_smoker = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_smoker.RData"))

GTEx_nonsmoker_genes = unlist(GTEx_nonsmoker$leadingEdge[GTEx_nonsmoker$pathway == pathway])
GTEx_smoker_genes = unlist(GTEx_smoker$leadingEdge[GTEx_smoker$pathway == pathway])

TCGA_nonsmoker_genes = unlist(TCGA_nonsmoker$leadingEdge[TCGA_nonsmoker$pathway == pathway])
TCGA_smoker_genes = unlist(TCGA_smoker$leadingEdge[TCGA_smoker$pathway == pathway])

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

# Ribosome pathway length
length(pathways[[pathway]])
ribosome_genes = pathways[[pathway]]

################################
# Check ribosomal gene expression: GTEx
# Load expression matrices
GTEx_TPM = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_TPM.txt"))
GTEx_TPM$V1 = gsub("\\..*","", GTEx_TPM$V1)
genes = GTEx_TPM$V1
GTEx_TPM = GTEx_TPM[,-1]
colnames(GTEx_TPM) = unlist(lapply(colnames(GTEx_TPM),function(x){substr(x,1,nchar(x)-2)}))
rownames(GTEx_TPM) = genes

GTEx_IDs = colnames(GTEx_TPM)
# Format GTEx IDs to match with phenotypic data: GTEx
GTEx_IDs = unlist(lapply(strsplit(GTEx_IDs, split=".", fixed = T),
                         function(x){paste(x[1:2], collapse= "-")}))

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(GTEx_TPM), gene_info.TCGA$gene_id)]

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

####### All ribosome genes #########
# Create dataframe for boxplot
# Get mean expression of ribosomal genes
ribosome_expression = colMeans(GTEx_TPM[which(gene_name %in% ribosome_genes),])

fulldata = data.frame(ribosome_expression)
fulldata$gender = factor(gender, levels = c("FEMALE", "MALE"))
fulldata$smoking_status = smoking_status

# Violinplot with jitters
maintitle = "GTEx: Distribution of Mean Expression of Ribosomal Genes"
ggplot(fulldata, aes(x = smoking_status, y = ribosome_expression, fill = gender)) + geom_violin() + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("Mean expression")
ggplot(fulldata, aes(x = smoking_status, y = ribosome_expression, fill = gender)) + geom_half_violin(side = "r") + geom_half_boxplot(side = "l") + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("Mean expression")

####### Ribosome genes exclusive to one group #########
# Get mean expression of genes exclusive to GTEx nonsmokers
exclusive_genes = setdiff(GTEx_nonsmoker_genes, GTEx_smoker_genes)
ribosome_expression = colMeans(GTEx_TPM[which(gene_name %in% exclusive_genes),])

# Create dataframe for boxplot
fulldata = data.frame(ribosome_expression)
fulldata$gender = factor(gender, levels = c("FEMALE", "MALE"))
fulldata$smoking_status = smoking_status

# Violinplot with jitters
maintitle = paste("GTEx: Mean expression distribution of ribosomal genes\nleading sex difference in nonsmokers but not in smokers") 
ggplot(fulldata, aes(x = smoking_status, y = ribosome_expression, fill = gender)) + geom_half_boxplot() + geom_half_violin(side = "r") + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("Mean expression")

####### Ribosome genes leading in one group #########
# Get mean expression of genes exclusive to GTEx nonsmokers
ribosome_expression = colMeans(GTEx_TPM[which(gene_name %in% GTEx_nonsmoker_genes),])

# Create dataframe for boxplot
fulldata = data.frame(ribosome_expression)
fulldata$gender = factor(gender, levels = c("FEMALE", "MALE"))
fulldata$smoking_status = smoking_status

# Violinplot with jitters
maintitle = paste("GTEx: Mean expression distribution of ribosomal genes\nleading sex difference in nonsmokers") 
ggplot(fulldata, aes(x = smoking_status, y = ribosome_expression, fill = gender)) + geom_half_boxplot() + geom_half_violin(side = "r") + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("Mean expression")

# Wilcoxon Test
nsf = fulldata$ribosome_expression[which(fulldata$smoking_status == "No" & fulldata$gender == "FEMALE")]
nsm = fulldata$ribosome_expression[which(fulldata$smoking_status == "No" & fulldata$gender == "MALE")]
sf = fulldata$ribosome_expression[which(fulldata$smoking_status == "Yes" & fulldata$gender == "FEMALE")]
sm = fulldata$ribosome_expression[which(fulldata$smoking_status == "Yes" & fulldata$gender == "MALE")]

wilcox.test(nsf, nsm, alternative = "greater")
wilcox.test(sf, sm, alternative = "less")

################################
# Check ribosomal gene expression: TCGA
# Load expression matrices
TCGA_TPM = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_TPM.txt"))
genes = gsub("\\..*","", TCGA_TPM$V1)
rownames(TCGA_TPM) = genes
TCGA_TPM = TCGA_TPM[,-1]

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(TCGA_TPM), gene_info.TCGA$gene_id)]

# remove X from the starting of TCGA IDs
colnames(TCGA_TPM)[which(substring(colnames(TCGA_TPM), 1,1) == "X")] = substring(colnames(TCGA_TPM)[which(substring(colnames(TCGA_TPM), 1,1) == "X")], 2, length(colnames(TCGA_TPM)[which(substring(colnames(TCGA_TPM), 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(colnames(TCGA_TPM), split=".", fixed = T),
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

####### All ribosome genes #########
# Get mean expression of ribosomal genes
ribosome_expression = colMeans(TCGA_TPM[which(gene_name %in% ribosome_genes),])

# Create dataframe for boxplot
fulldata = data.frame(ribosome_expression)
fulldata$gender = factor(gender, levels = c("FEMALE", "MALE"))
fulldata$smoking_status = smoking_status

# Violinplot with jitters
maintitle = "TCGA: Distribution of Mean Expression of Ribosomal Genes"
ggplot(fulldata, aes(x = smoking_status, y = ribosome_expression, fill = gender)) + geom_violin() + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("Mean expression")
ggplot(fulldata, aes(x = smoking_status, y = ribosome_expression, fill = gender)) + geom_half_violin(side="r") + geom_half_boxplot() + geom_jitter(size=0.25) + theme_bw() + ggtitle(maintitle) + xlab("Smoking Status") + ylab("Mean expression")

########################################
# Violinplot for Ribosomal Gene Targeting in GTEx
nonsmoker_limma = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_GTEx_nonsmoker.txt"))
gene_name = gene_info.TCGA$gene_name[match(nonsmoker_limma$V1, gene_info.TCGA$gene_id)]
ribosome_nonsmoker = nonsmoker_limma[which(gene_name %in% ribosome_genes),-1]

smoker_limma = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_GTEx_smoker.txt"))
gene_name = gene_info.TCGA$gene_name[match(smoker_limma$V1, gene_info.TCGA$gene_id)]
ribosome_smoker = smoker_limma[which(gene_name %in% ribosome_genes),-1]

library(ggplot2)
library(ggrepel)

tb = ribosome_nonsmoker
# cutoff for logFC
FClim = 100
# add a column of NAs
tb$direction <- "Not significant"
# if log2Foldchange > 25 and -log10(pvalue) > 5, set as "UP" 
tb$direction[tb$logFC > FClim & tb$P.Value < 0.05] <- "Female"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tb$direction[tb$logFC < -FClim & tb$P.Value < 0.05] <- "Male"

mycolors <- c("red", "blue", "grey")
names(mycolors) <- c("Female", "Male", "Not significant")

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
tb$delabel <- NA
tb$delabel[tb$direction != "Not significant"] <- tb$gene_name[tb$direction != "Not significant"]

# Remove gene names that copntain "-" or "." in name
tb$delabel[grep("-",tb$delabel)] <- NA
tb$delabel[grep(".",tb$delabel, fixed = T)] <- NA

# Re-plot but this time color the points with "diffexpressed"
ggplot(data=tb, aes(x=logFC, y=-log10(P.Value), col=direction, label = delabel)) +
  geom_point() + theme_minimal() +
  geom_vline(xintercept=c(-FClim, FClim), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_colour_manual(values = mycolors) +
  ggtitle("GTEx Nonsmoker: Volcanoplot of Differentially Targeted Ribosomal Genes by Sex") +
  geom_text_repel(size = 3)


##################################
# Ribosomal gene heatmap of t-statistics
pathway = "KEGG_RIBOSOME"

GTEx_nonsmoker = read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_GTEx_nonsmoker.txt", sep="")
GTEx_smoker = read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_GTEx_smoker.txt", sep="")

TCGA_nonsmoker = read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_TCGA_nonsmoker.txt", sep="")
TCGA_smoker = read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/limma_TCGA_smoker.txt", sep="")

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

# Ribosome pathway length
length(pathways[[pathway]])
ribosome_genes = pathways[[pathway]]

ribosoma_ids = gsub("\\..*","",gene_info.TCGA$gene_id[match(ribosome_genes, gene_info.TCGA$gene_name)])

tab = cbind(GTEx_nonsmoker$t[match(ribosoma_ids, rownames(GTEx_nonsmoker))],
            GTEx_smoker$t[match(ribosoma_ids, rownames(GTEx_smoker))],
            TCGA_nonsmoker$t[match(ribosoma_ids, rownames(TCGA_nonsmoker))],
            TCGA_smoker$t[match(ribosoma_ids, rownames(TCGA_smoker))])

colnames(tab) = c("GTEx_nonsmoker", "GTEx_smoker", "TCGA_nonsmoker", "TCGA_smoker")
rownames(tab) = ribosome_genes

tab = tab[complete.cases(tab),]

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/heatmaps/Ribosomal_genes_tStat.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="tStat", Colv = F)
dev.off()


