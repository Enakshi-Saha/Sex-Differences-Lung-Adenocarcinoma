setwd("~/OneDrive/Documents/Channing/Projects/cancer/lung_cancer/ALA/GDSC_codeReview")

library(Biobase)
library(gplots)
library(ggfortify)
eset <- get(load("dataProcessed/gdsc_exp_rma.RData"))
# Subset the ExpressionSet for Y chromosome genes
k <- which(fData(eset)$chr == "Y")
esetY <- eset[k,]
esetY
exp <- exprs(esetY)
SD <- apply(exp,1,sd)
names(SD)[SD>2]
# Select only Y genes with SD > 2 (7 genes)
exp2 <- exp[names(SD)[SD>2],]
hist(colSums(exp2))
hist(colMeans(exp2))

heatmap.2(exp2,trace='none')
myPCA <- prcomp(t(exp2))
autoplot(myPCA, scale=FALSE, data=pData(esetY),  main="Y chromosome genes expression")

# Annotate sex based on Y  chr genes mean expression
sex <- rep("remove",ncol(exp2))
sex[colMeans(exp2) > 5] <- "male"
sex[colMeans(exp2) < 3.8] <- "female"
table(sex)
pData(esetY)$sex <- sex
myPCA <- prcomp(t(exp2))
autoplot(myPCA, scale=FALSE,data=pData(esetY), colour="sex",  main="Y chromosome genes expression")

k <- which(pData(esetY)$sex == "remove")
objplot <- esetY[rownames(exp2),-k]
cols <- pData(objplot)$sex
heatmapColColors <- c("red", "blue")[factor(cols)]
heatmap.2(exprs(objplot),trace='none',labRow=fData(objplot)$hgnc, labCol="",ColSideColors=heatmapColColors,margins=c(10,12))
legend(3.5,4, legend=unique(cols),fill=unique(heatmapColColors), xpd=TRUE, box.lwd=NA, cex=0.7) 


# Save sex annotation to eset
pData(eset)$sex <- sex
table(pData(eset)$sex)

colnames(pData(eset))[which(colnames(pData(eset))=="Factor.Value.cell_line.")] <- "cell_line"

save(eset,file="dataProcessed/gdsc_exp_rma_sex.RData")
