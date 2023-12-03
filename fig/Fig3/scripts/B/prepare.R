library(data.table)
library(dplyr)
library(stats)
library(ggpubr)
library(caret)
library(precrec)
library(readxl)
library(EnvStats)
library(gghalves)
library(ggdist)
library(reshape2)

directory = "/path/to/data"

input_dir=paste0(directory, "/data/dataInput.old")
output_dir=paste0(directory, "/fig/Fig2/data/B")

dt=as.data.frame(fread(paste0(input_dir, "/bigTable.mergedAll.tsv")))
row.names(dt)=dt$probeID

set1 <- read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/data/train/DMP/TB-HC/BMA.txt"), header=T)$probe_ID[-1]
set2 <- read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/data/train/DMP/TB-HC/RF_min.txt"), header=T)$name
set3 <- read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/data/train/DMP/TB-HC/PLS_min.txt"), header=T)$name
set4 <- read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/data/train/DMP/TB-HC/LASSO.txt"), header=T)$probeID

SigCpG=union(union(union(set1, set2), set3), set4)

Sample=read.table(paste0(directory, "/config/design_train_data.tsv"),sep="\t", header=T)
names(Sample)=c("Accession", "Class")

df=dt[dt$probeID%in%SigCpG, Sample$Accession]
tdf=data.frame(t(df))
tdf$sample_ID=rownames(tdf)
tdf$Class=merge(tdf, Sample, by.x="sample_ID", by.y="Accession")$Class
trainDT=arrange(tdf, desc(Class))
trainDT=trainDT[,-grep("sample_ID", names(trainDT))]

train_plot=melt(trainDT)

write.table(train_plot, paste0(output_dir, "/boxplotData.tsv"), quote=F, row.names=F, sep="\t")

### prepare limma p-value
dt=as.data.frame(fread(paste0(directory, "/data/train/DMP.current/TB-HC/DMPs.hg38.tsv")))
dt=dt[dt$probeID%in%SigCpG,c("probeID", "ind.fdr")]
write.table(dt, paste0(output_dir, "/Limma_indFDR.tsv"), quote=F, row.names=F, sep="\t")

