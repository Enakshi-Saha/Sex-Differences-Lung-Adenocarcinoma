library("data.table")

# Get Y genes
# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
Y_genes = gene_info.TCGA$gene_name[which(gene_info.TCGA$seqnames == "chrY")]

expression = fread("/home/esaha/validation_data/GSE68465/GSE68465_expression.txt")
expression = data.frame(expression)

# Transform log2(1+value)
expression[, -c(1,2)] = log2(1+expression[, -c(1,2)])

phenotypes = data.frame(fread("/home/esaha/validation_data/GSE68465/GSE68465_phenotype.txt"))

# Check sex-misannotation
gender = phenotypes$Sex.ch1[match(colnames(expression)[-(1:2)], phenotypes$geo_accession)] 
library(ggplot2)
#### Remove "//" from gene_symbol and only keep the first name
gene_symbols = expression[,2]
gene_symbols2 = lapply(strsplit(gene_symbols, split = "//", fixed= T), function(x){gsub(" ", "", x[1])})
Y_genes_subsample = gene_symbols[which(unlist(gene_symbols2) %in% Y_genes)]

# PCA plot to detect sex mis-specification
Y_expr = expression[which(gene_symbols %in% Y_genes_subsample),-c(1,2,which(is.na(gender))+2)]
gender = factor(as.character(gender[-which(is.na(gender))]))
pca = prcomp(t(Y_expr))
U <- pca$x
theme_set(bigstatsr::theme_bigstatsr(0.8))
qplot(U[, 1], U[, 2], main = "PCA of Y genes: GSE68465", colour = gender) + coord_equal() + xlab("PC1") +ylab("PC2")

# Samples to exclude because of sex-misannotation
# U[,1] refers to the first principal component of the expression of all genes on the Y chromosome. 
# From a scatterplot of the first two principal components we observed 
# that for almost all female samples the values of the first principal component were lager than 0 
# and for almost all male samples the values were smaller than 2.5. 
# The samples that violated these conditions were excluded from the analysis as they clustered with the opposite sex. 
female_wrong = colnames(Y_expr)[which(U[,1] <0 & gender == "Female")]
male_wrong = colnames(Y_expr)[which(U[,1] > 2.5 & gender == "Male")]
sex_unknown = phenotypes$V1[which(is.na(phenotypes$Sex.ch1))]

expression = expression[,-which(colnames(expression) %in% c(female_wrong, male_wrong, sex_unknown))]
phenotypes = phenotypes[match(colnames(expression)[-c(1,2)], phenotypes$V1),]

sample_IDs = colnames(expression)[-c(1,2)]
sample_source = phenotypes$source_name_ch1[match(sample_IDs, phenotypes$geo_accession)]

pca = prcomp(t(data.matrix(expression[,-c(1,2)])))
U <- pca$x
plot(U[, 1], U[, 2], main = "PCA of GSE68465", col = factor(sample_source), pch = 19, xlab = "PC1", ylab = "PC2")

# Batch Correction
library(sva)
expression_corrected <- ComBat(expression[,-c(1,2)], batch = sample_source)
rownames(expression_corrected) = expression[,1]

pca_corrected = prcomp(t(data.matrix(expression_corrected)))
U <- pca_corrected$x
plot(U[, 1], U[, 2], main = "PCA of GSE68465 (batch-corrected)", col = factor(sample_source), pch = 19, xlab = "PC1", ylab = "PC2")

# For every gene symbol, compute the standard deviation (SD) for all the probes. 
# Keep the probe with highest SD and discard the rest. This gives us gene expression matrix for 13516 genes.
probewise_sd = unlist(apply(expression_corrected, MARGIN = 1, sd))
df = data.frame(cbind(gene_symbols, probewise_sd, names(probewise_sd) ))
# df %>% group_by(gene_symbols) %>% top_n(1, probewise_sd)
library(dplyr)
filtered = df %>% group_by(gene_symbols) %>% dplyr::slice(which.max(probewise_sd))
genes_filtered = filtered$gene_symbols
probes_filtered = filtered$V3
expression_filtered = expression_corrected[which(rownames(expression_corrected) %in% probes_filtered),]
rownames(expression_filtered) = genes_filtered[match(rownames(expression_filtered), probes_filtered)]

# Set female Y gene expression to "NA"
female_id = phenotypes$geo_accession[which(phenotypes$Sex.ch1 == "Female")]
male_id = phenotypes$geo_accession[which(phenotypes$Sex.ch1 == "Male")]

expression_male = expression_filtered[,which(colnames(expression_filtered) %in% male_id)]
expression_female = expression_filtered[,which(colnames(expression_filtered) %in% female_id)]

expression_female[which(rownames(expression_female) %in% Y_genes_subsample),] = "NA"


write.table(expression_male, file = "/home/esaha/validation_data/validation_newdata2023/GSE68465_male.txt", 
            row.names = T, col.names = T)
write.table(expression_female, file = "/home/esaha/validation_data/validation_newdata2023/GSE68465_female.txt", 
            row.names = T, col.names = T)

write.table(expression_filtered, file = "/home/esaha/validation_data/validation_newdata2023/GSE68465_expression_filtered.txt", 
            row.names = T, col.names = T)
write.table(phenotypes, file = "/home/esaha/validation_data/validation_newdata2023/GSE68465_phenotypes.txt", 
            row.names = T, col.names = T)

