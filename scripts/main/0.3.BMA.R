library(data.table)
library(dplyr)
library(BMA)
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
x <- t.bmatrix[,sig.cpg]
y <- ss$Class

BMAmodel <- bicreg(x, y, OR=20, strict = TRUE)
summary(BMAmodel)
imageplot.bma(BMAmodel)

B <- data.frame(coeff=BMAmodel$ols[1,])
B$probe_ID <- rownames(B)
B2 <- B[!(B$coeff==0),]

### wilcoxon ranksum
t.bmatrix$sN=rownames(t.bmatrix)
for (i in B2[grep("cg", B2$probe_ID), "probe_ID"]) {
tmp=merge(t.bmatrix[,c(i, "sN")], ss, by.x="sN", by.y="Accession")

stat=wilcox.test(x=tmp[tmp$Class==0,i], y=tmp[tmp$Class==1,i])
B2[B2$probe_ID==i,"wilcox"]=stat$p.value
}


write.table(B2, file = paste0(data_dir,"/data/train/DMP/TB-HC/BMA.txt"), sep="\t", quote=F, row.names=F)
