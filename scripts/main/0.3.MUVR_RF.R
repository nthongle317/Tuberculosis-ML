# installation
# install.packages('remotes')
# library(remotes)
# install_gitlab('CarlBrunius/MUVR')

# load
library(data.table)
library(dplyr)
library(ggplot2)
library(doParallel)
library(MUVR)

set.seed(123)

data_dir="/path/to/data_dir"

### Get phenotype from sample sheet
ss=as.data.frame(fread(paste0(data_dir, "/config/train.tsv")))
table(ss$status)
sN=ss$Accession

### Get probe beta value of each patient
dt=as.data.frame(fread(paste0(data_dir, "/data/dataInput.old/bigTable.mergedAll.tsv")))
cpg=as.data.frame(fread(paste0(data_dir,"/data/train/DMP/TB-HC/DMPs.hg38.tsv")))
sig.cpg=cpg[abs(cpg$diff)>=0.2,"probeID"]

bmatrix=dt[dt$probeID%in%sig.cpg,c(sN, "probeID")]
list.cpg=bmatrix$probeID
rownames(bmatrix)=list.cpg
t.bmatrix=as.data.frame(t(bmatrix[sig.cpg,-grep("probeID", names(bmatrix))]))

### Create input for MUVR
x=as.matrix(t.bmatrix[,sig.cpg])
# Response
y=ss$Class
# Subject identifier
z=ss$Accession

### RF (100 repetition)
nCore=5
# Number of MUVR repetitions - differ nRep -> differ results 
nRep=100
# Number of outer cross-validation segments          
nOuter= 8   
# Proportion of variables kept per iteration            
varRatio=0.8    

#####################RF
# Selected core modelling algorithm        
method='RF'     
#MUVR
cl=makeCluster(nCore) 
registerDoParallel(cl)
classModel=MUVR(X=x, Y=y, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)
stopCluster(cl)

saveRDS(classModel, file=paste0(data_dir, "/data/train/DMP/TB-HC/MUVR_RF_result.RDS"))

### Choose model (min, mid, max)
pdf(file=paste0(data_dir,"/data/train/DMP/TB-HC/MUVR_RF_plot.pdf"), width=20, height=10)
plotVIP(classModel, model='min')
plotVIP(classModel, model='mid')
plotVIP(classModel, model='max')
dev.off()

# write table
MUVR_res_min = getVIP(classModel, model='min')
# MUVR_res_min
write.table(MUVR_res_min, paste0(data_dir, "/data/train/DMP/TB-HC/RF_min.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

MUVR_res_mid = getVIP(classModel, model='mid')
# MUVR_res_mid
write.table(MUVR_res_mid, paste0(data_dir, "/data/train/DMP/TB-HC/RF_mid.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

MUVR_res_max = getVIP(classModel, model='max')
# MUVR_res_max
write.table(MUVR_res_max, paste0(data_dir, "/data/train/DMP/TB-HC/RF_max.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
