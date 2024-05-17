# Make boxplot of ln(IC50) for candidate drugs by sex

library(Biobase)
library(gghalves)
library(ggplot2)

setwd("/home/ubuntu/GDSC_codeReview")
tab1 <- read.delim("dataDownload/GDSC1_fitted_dose_response_25Feb20.txt",stringsAsFactors = F)
tab2 <- read.delim("dataDownload/GDSC2_fitted_dose_response_25Feb20.txt",stringsAsFactors = F)
tab <- rbind(tab1,tab2)

###
#GDSC data
eset <- get(load("dataProcessed/gdsc_exp_rma_sex.RData"))
annot <- pData(eset)[match(tab$CELL_LINE_NAME,pData(eset)$cell_line),]

table(pData(eset)$cell_line %in% tab$CELL_LINE_NAME)
table(tab$CELL_LINE_NAME %in% pData(eset)$cell_line)

# there are some cell lines the name does not match between expression and sensitivity dta. So to match them, I removed "-" and made lower caps
tab$cell_line_mod <-  tolower(tab$CELL_LINE_NAME)
tab$cell_line_mod <- sub("786-0","786-o",tab$cell_line_mod)
tab$cell_line_mod <- sub("-","",tab$cell_line_mod)
tab$cell_line_mod <- sub("\\.","",tab$cell_line_mod)

pData(eset)$cell_line_mod <-  tolower(pData(eset)$cell_line)
pData(eset)$cell_line_mod <- sub("-","",pData(eset)$cell_line_mod)
pData(eset)$cell_line_mod <- sub("\\.","",pData(eset)$cell_line_mod)

sort(unique(pData(eset)$cell_line_mod))

table(unique(tab$cell_line_mod) %in% pData(eset)$cell_line_mod)
# FALSE  TRUE 
# 49   937 

annot <- pData(eset)[match(tab$cell_line_mod,pData(eset)$cell_line_mod),]
tab <- cbind(tab,"sex_Y"=annot[,"sex"])

# Cell model passaports from sanger
cell_passp <- read.csv("dataDownload/model_list_20220315.csv")
annot2 <- cell_passp[match(tab$CELL_LINE_NAME,cell_passp$model_name),]
annot2$sex_origin <- tolower(annot2$gender)
tab <- cbind(tab,"sex_origin"=annot2$sex_origin)

# create variable that matches sex at origin and Y expression
tab$sex_matchYorigin <- ifelse(tab$sex_Y==tab$sex_origin, yes=tab$sex_Y, no="inconsistent")


alldrug <- c("Trametinib", "Panobinostat", "Vorinostat", "Dactinomycin")

pdf("/home/esaha/paper_plots/violinplot/paper_LUAD_violinplot_alldrugs_lnIC50_noOutlier.pdf",h=5,w=3)
for(drug in alldrug){
  tab2 <- tab[which(tab$DRUG_NAME == drug),]
  # remove sex-sepecific cell lines
  cl <- c("PRAD", "TGCT", "UCEC", "OV", "CESC", "UCS", "BRCA")
  cl <- which(tab2$TCGA_DESC %in% cl)
  tab2 <- tab2[-cl,]
  
  
  # Combine same cell lines (technical replicates) by median of IC50
  sens <- by(tab2,tab2$CELL_LINE_NAME, function(x) median(x[,"LN_IC50"]))
  sens <- as.matrix(sens)
  colnames(sens) <- "LN_IC50"
  sens <- as.data.frame(sens)
  
  sex_info <- tab2[match(rownames(sens), tab2$CELL_LINE_NAME),c("sex_Y", "sex_origin", "sex_matchYorigin")]
  sens <- cbind(sens,sex_info)
  tab2 <- sens
  
  k <- which(tab2$sex_matchYorigin %in% c("male", "female"))
  tab2 <- tab2[k,]
  
  res <- wilcox.test(tab2$LN_IC50 ~ tab2$sex_matchYorigin)
  
  Title2 <- paste0(drug,"\n","p=",round(res$p.value,5))
  #boxplot(tab2$LN_IC50 ~ tab2$sex_matchYorigin,ylab="ln(IC50)",main=Title2,col=c("red","blue"),xlab=NA,outline=F)
  
  # violinplot
  tab2$sex = tab2$sex_matchYorigin
  violinplot = ggplot(tab2, aes(x = sex_matchYorigin, y = LN_IC50, fill = sex)) + geom_violin() + geom_jitter(show.legend = FALSE) + theme_bw() + ggtitle(Title2) + xlab("") + ylab("ln(IC50)") + theme(legend.position = "none")
  print(violinplot)
  halfviolinplot = ggplot(tab2, aes(x = sex_matchYorigin, y = LN_IC50, fill = sex)) + geom_half_violin(side="r") + geom_half_boxplot() + geom_jitter(show.legend = FALSE, size = 0.5) + theme_bw() + ggtitle(Title2) + xlab("") + ylab("ln(IC50)") + theme(legend.position = "none")
  print(halfviolinplot)
}
plot.new()
ggplot(tab2, aes(x = sex_matchYorigin, y = LN_IC50, fill = sex)) + geom_violin() + geom_jitter(show.legend = FALSE) + theme_bw() + ggtitle(Title2) + xlab("") + ylab("ln(IC50)")
# legend("top",c("Female","Male"),fill=c("red","blue"),box.lty = 0)
dev.off()


# Find the number of total cell lines used in the analysis after removing sex-specific tumors
cl <- c("PRAD", "TGCT", "UCEC", "OV", "CESC", "UCS", "BRCA")
cl <- which(tab$TCGA_DESC %in% cl)
tab3 <- tab[-cl,]
k <- which(tab3$sex_matchYorigin %in% c("male", "female"))
tab3 <- tab3[k,]
# remove cell line duplicated entry
tab3 <- tab3[!(duplicated(tab3$CELL_LINE_NAME)),]
table(tab3$sex_matchYorigin)
# female   male 
# 264    227 
wilcox.test(tab3$LN_IC50 ~ tab3$sex_matchYorigin)

