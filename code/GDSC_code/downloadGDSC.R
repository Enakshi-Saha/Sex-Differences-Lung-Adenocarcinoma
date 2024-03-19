# Download GDSC expression data from ArrayExpress (E-MTAB-3610), RMA processing, annotate probe ID and chromosome location
# A-GEOD-13667 - [HG-U219] Affymetrix Human Genome U219 Array
# New eset object: "dataProcessed/gdsc_exp_rma.RData"

setwd("~/OneDrive/Documents/Channing/Projects/cancer/lung_cancer/ALA/GDSC_codeReview")

library(ArrayExpress)
library(biomaRt)
library(oligo)

setwd("/Users/camila/OneDrive/Documents/Channing/Projects/cancer/lung_cancer/ALA/GDSC_codeReview")

# Dowload and save raw expression 
rawset = ArrayExpress("E-MTAB-3610")
save(rawset,file="dataDownload/E-MTAB-3610.RData")

# Load the raw expression saved and normalize by RMA
obj <- get(load("dataDownload/E-MTAB-3610.RData"))
obj2 <- oligo::rma(obj)

# There are duplicated cell lines, so cannot use them as colnames for expression
pData(obj2)$Characteristics.cell.line.[duplicated(pData(obj2)$Characteristics.cell.line.)]
#[1] "PC-3"      "SK-MEL-28" "KM-H2"     "OACp4C"    "OCI-AML5" 

# To extract the processed expression matrix:
exp <- exprs(obj2)[1:nrow(exprs(obj2)),1:ncol(exprs(obj2))]
# convert probe id to gene symbol (I downloaded the platform information from ArrayExpress)
annot <- read.delim("dataDownload/A-GEOD-13667.adf.txt",stringsAsFactors = F, skip=15)
genes <- rownames(exp)
genbank <- annot[match(genes,annot[,1]),2]
# Keep only the first refseq ID for the  probes that have more than one matching gene 
genbank <- unlist(lapply(strsplit(genbank,";"),function(x) x[1]))
# Add symbol ID, ENS, chr
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",host = "www.ensembl.org")
annot2 <- getBM(attributes=c('refseq_mrna', 'hgnc_symbol', 'ensembl_gene_id', 'chromosome_name'), 
               filters = 'refseq_mrna', 
               values = genbank, 
               mart = ensembl)
ftab <- annot2[match(genbank, annot2$refseq_mrna),]
colnames(ftab) <- c("refseq", "hgnc", "ens", "chr")
rownames(ftab) <- genes

# Make eset object
eset <- ExpressionSet(assayData=as.matrix(exp), 
                 phenoData=AnnotatedDataFrame(pData(obj2)),
                 featureData=AnnotatedDataFrame(ftab))
save(eset,file="dataProcessed/gdsc_exp_rma.RData")



