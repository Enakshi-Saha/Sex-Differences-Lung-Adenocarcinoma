##### Read lioness networks:
rm(list=ls())
dataset = "TCGA"
gender = "female"

library(data.table)
setwd(paste("/home/ubuntu/", dataset, "_Lung_lioness_", gender, "/", sep = ""))
files_ls = list.files(path=".", pattern=".txt", all.files=TRUE,
                      full.names=TRUE)

indices = unlist(lapply(files_ls, function(x){unlist(strsplit(x, split = ".", fixed = T))[3]}))

# Read TF_gene names from PANDA networks
panda = read.table(paste("/home/ubuntu/", dataset, "_lung_panda_", gender, ".txt", sep = ""), quote="\"", comment.char="")
panda = data.frame(panda)

indegree = aggregate(panda$V4, by=list(Category=panda$V2), FUN=sum)

genes = indegree[,1]

data = list()
i=0
for (file_name in files_ls){
  i = i+1
  new_data = fread(file_name)
  data[[i]] = colSums(new_data)
  cat(paste(i,"\n"))
}

lioness = do.call(cbind, data)

rownames(lioness) = genes
colnames(lioness) = indices

head(lioness)
dim(lioness)

setwd("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/")
write.table(lioness, file=paste(dataset, "_lung_lioness_indegree_", gender, ".txt", sep = ""), row.names=T, col.names=T)

##############################################

##### Read lioness networks:
rm(list=ls())
dataset = "GSE68465"
gender = "female"

library(data.table)
setwd(paste("/home/ubuntu/", "validation_lioness_", gender, "_", dataset, "/", sep = ""))
files_ls = list.files(path=".", pattern=".txt", all.files=TRUE,
                      full.names=TRUE)

indices = unlist(lapply(files_ls, function(x){unlist(strsplit(x, split = ".", fixed = T))[3]}))

# Read TF_gene names from PANDA networks
panda = read.table(paste("/home/ubuntu/", "validation_panda_", gender, "_", dataset, ".txt", sep = ""), quote="\"", comment.char="")
panda = data.frame(panda)

indegree = aggregate(panda$V4, by=list(Category=panda$V2), FUN=sum)

genes = indegree[,1]

data = list()
i=0
for (file_name in files_ls){
  i = i+1
  new_data = fread(file_name)
  data[[i]] = colSums(new_data)
  cat(paste(i,"\n"))
}

lioness = do.call(cbind, data)

rownames(lioness) = genes
colnames(lioness) = indices

head(lioness)
dim(lioness)

setwd("/home/esaha/validation_data/validation_newdata2023/validation_indegree/")
write.table(lioness, file=paste(dataset, "_lioness_indegree_", gender, ".txt", sep = ""), row.names=T, col.names=T)
