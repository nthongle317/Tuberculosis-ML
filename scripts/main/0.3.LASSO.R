library(data.table)
library(dplyr)
library(glmnet)
set.seed(123)

data_dir="/path/to/data_dir"

### Get phenotype from sample sheet
ss=as.data.frame(fread(paste0(data_dir, "/config/train.tsv")))
ss$Class[ss$Class == "TB"] <- 1
ss$Class[ss$Class == "HC"] <- 0

table(ss$Class)
sN=ss$Accession

### Get probe beta value of each patient
dt=as.data.frame(fread(paste0(data_dir, "/data/dataInput/bigTable.mergedAll.tsv")))
cpg=as.data.frame(fread(paste0(data_dir,"/data/train/DMP/TB-HC/DMPs.hg38.tsv")))
sig.cpg=cpg[abs(cpg$diff)>=0.2,"probeID"]

bmatrix=dt[dt$probeID%in%sig.cpg,c(sN, "probeID")]
list.cpg=bmatrix$probeID
rownames(bmatrix)=list.cpg
t.bmatrix=as.data.frame(t(bmatrix[sig.cpg,-grep("probeID", names(bmatrix))]))
### Prepare for BMA model selection
x <- as.matrix(t.bmatrix[,sig.cpg])
y <- as.numeric(ss$Class)

cvfit <- cv.glmnet(x, y, family="gaussian", alpha=1, nlambda=100)
plot(cvfit)
saveRDS(cvfit, file = paste0(data_dir,"HieuPackage/LASSO_result.rds"))

A  <- coef(cvfit, s = "lambda.1se")
B <- as.data.frame(as.matrix(A))
B$probeID <- rownames(B)
B <- B[-1,]
colnames(B)[1] <- "coeff" 
B2 <- B[!(B$coeff==0),]

write.table(B2, file = paste0(data_dir,"/data/train/DMP/TB-HC/LASSO.txt"), sep="\t", quote=F, row.names=F)