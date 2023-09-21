#### Heatmap 1: GTEx Smoker vs Nonsmoker ####

nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_nonsmoker.RData"))
smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_smoker.RData"))

smoker = smoker[match(nonsmoker$pathway, smoker$pathway),]

tab = cbind(nonsmoker$NES, smoker$NES)

rownames(tab) = nonsmoker$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((nonsmoker$padj<0.05) | (smoker$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("nonsmoker", "smoker")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/heatmaps/GTEx_smoker_nonsmoker_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

#### Heatmap 2: GTEx vs TCGA Nonsmoker ####

gtex_nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_nonsmoker.RData"))
gtex_smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_smoker.RData"))
tcga_nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_nonsmoker.RData"))
tcga_smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_smoker.RData"))

gtex_smoker = gtex_smoker[match(gtex_nonsmoker$pathway, gtex_smoker$pathway),]
tcga_nonsmoker = tcga_nonsmoker[match(gtex_nonsmoker$pathway, tcga_nonsmoker$pathway),]
tcga_smoker = tcga_smoker[match(gtex_nonsmoker$pathway, tcga_smoker$pathway),]

tab = cbind(gtex_nonsmoker$NES, tcga_nonsmoker$NES)

rownames(tab) = gtex_nonsmoker$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((tcga_nonsmoker$padj<0.05) | (tcga_smoker$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("GTEx", "TCGA")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/heatmaps/GTEx_TCGA_nonsmoker_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

#### Heatmap 3: GTEx vs TCGA Smoker ####

gtex_nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_nonsmoker.RData"))
gtex_smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_smoker.RData"))
tcga_nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_nonsmoker.RData"))
tcga_smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_smoker.RData"))

gtex_smoker = gtex_smoker[match(gtex_nonsmoker$pathway, gtex_smoker$pathway),]
tcga_nonsmoker = tcga_nonsmoker[match(gtex_nonsmoker$pathway, tcga_nonsmoker$pathway),]
tcga_smoker = tcga_smoker[match(gtex_nonsmoker$pathway, tcga_smoker$pathway),]

tab = cbind(gtex_smoker$NES, tcga_smoker$NES)

rownames(tab) = gtex_nonsmoker$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((tcga_nonsmoker$padj<0.05) | (tcga_smoker$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("GTEx", "TCGA")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/heatmaps/GTEx_TCGA_smoker_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

###### Heatmap 4: GSE68465 Smoker-Nonsmoker ####

tcga_nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_nonsmoker.RData"))
tcga_smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_smoker.RData"))

tcga_smoker = tcga_smoker[match(tcga_nonsmoker$pathway, tcga_smoker$pathway),]

GSE68465_nonsmoker = get(load("~/validation_data/validation_newdata2023/validation_limma_GSEA/gsea_GSE68465_nonsmoker.RData"))
GSE68465_smoker = get(load("~/validation_data/validation_newdata2023/validation_limma_GSEA/gsea_GSE68465_smoker.RData"))

GSE68465_nonsmoker = GSE68465_nonsmoker[match(tcga_nonsmoker$pathway, GSE68465_nonsmoker$pathway),]
GSE68465_smoker = GSE68465_smoker[match(tcga_nonsmoker$pathway, GSE68465_smoker$pathway),]

tab = cbind(GSE68465_nonsmoker$NES, GSE68465_smoker$NES)

rownames(tab) = GSE68465_nonsmoker$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((tcga_nonsmoker$padj<0.05) | (tcga_smoker$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("Nonsmoker", "Smoker")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/heatmaps/GSE68465_smoker_nonsmoker_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()


###### Heatmap 5: GSE68465 Stages ####

tcga_nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_nonsmoker.RData"))
tcga_smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_TCGA_smoker.RData"))

tcga_smoker = tcga_smoker[match(tcga_nonsmoker$pathway, tcga_smoker$pathway),]

GSE68465_stageI = get(load("~/validation_data/validation_newdata2023/validation_limma_GSEA/gsea_GSE68465_stageI.RData"))
GSE68465_stageII = get(load("~/validation_data/validation_newdata2023/validation_limma_GSEA/gsea_GSE68465_stageII.RData"))
GSE68465_stageIII = get(load("~/validation_data/validation_newdata2023/validation_limma_GSEA/gsea_GSE68465_stageIII.RData"))

GSE68465_stageI = GSE68465_stageI[match(tcga_nonsmoker$pathway, GSE68465_stageI$pathway),]
GSE68465_stageII = GSE68465_stageII[match(tcga_nonsmoker$pathway, GSE68465_stageII$pathway),]
GSE68465_stageIII = GSE68465_stageIII[match(tcga_nonsmoker$pathway, GSE68465_stageIII$pathway),]

tab = cbind(GSE68465_stageI$NES, GSE68465_stageII$NES, GSE68465_stageIII$NES)

rownames(tab) = GSE68465_stageI$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((tcga_nonsmoker$padj<0.05) | (tcga_smoker$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("StageI", "StageII", "StageIII-IV")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/heatmaps/GSE68465_stages_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

######################################
#### Heatmap 5: LTRC Smoker vs Nonsmoker ####

GTEx_nonsmoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_nonsmoker.RData"))
GTEx_smoker = get(load("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/limma_GSEA/gsea_GTEx_smoker.RData"))

GTEx_smoker = GTEx_smoker[match(GTEx_nonsmoker$pathway, GTEx_smoker$pathway),]

GTEx_LTRC_NES_allPathways <- read.csv("~/paper_plots/heatmaps/GTEx_LTRC_smoker_heatmap_NES_allPathways.txt", sep="")

LTRC_NES_allPathways = GTEx_LTRC_NES_allPathways[match(GTEx_nonsmoker$pathway, rownames(GTEx_LTRC_NES_allPathways)),c(2,4)]

tab = LTRC_NES_allPathways

rownames(tab) = GTEx_nonsmoker$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((GTEx_nonsmoker$padj<0.05) | (GTEx_smoker$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("nonsmoker", "smoker")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

tab_subset = na.omit(tab_subset)

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/paper_plots/plots2023/heatmaps/LTRC_smoker_nonsmoker_heatmap.pdf",h=25,w=16)
heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

