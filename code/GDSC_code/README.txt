Code for Figure 8. 
Using GDSC drug sensitivity data to validate CLUEreg candidate genes.

Author: Camila M. Lopes-Ramos


1. downloadGDSC.R
This script downloads GDSC expression data from ArrayExpress (E-MTAB-3610) and normalize the data by RMA
Output: "dataProcessed/gdsc_exp_rma.RData"

2. sex_annotation.R
This script reads "dataProcessed/gdsc_exp_rma.RData" and add a new variable "sex", which uses Y  chr genes mean expression to define males and females.
Output: "dataProcessed/gdsc_exp_rma_sex.RData"

3. paper_LUAD_drugResponse.R
This script reads the GDSC expression set "dataProcessed/gdsc_exp_rma_sex.RData" and GDSC sensitivity data downloaded on 25Feb20 and stored on "dataDownload".
- I also considered the cell line sex of the individual origin. GDSC website points to sanger cell passport
model_list_20220315.csv
I downloaded this file from https://cellmodelpassports.sanger.ac.uk/downloads
- I only kept cell lines that had a match between sex from origin and sex defined by Y chr expression. This removed inconsistencies and cell lines that lost Y chromosome. 
- Make boxplot of ln(IC50) for candidate drugs by sex
output: "figure/paper_LUAD_boxplot_alldrugs_lnIC50_noOutlier.pdf" 
