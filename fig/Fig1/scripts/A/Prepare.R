library(data.table)
library(dplyr)

### Setup directory
directory = "/path/to/data"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/fig/Fig1/data/A")

### get beta matrix
dt <- as.data.frame(fread(paste0(input_dir, "/bigTable.mergedAll.tsv")))
row.names(dt)=dt$probeID
sampleSheet=read.table(paste0(directory, "/config/train.tsv"),sep="\t", header=T)
cpgList=read.table(paste0(directory, "/data/train/DMP/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)

SigCpG=cpgList[,"probeID"]

cpgList=cpgList[cpgList$probeID%in%SigCpG, ]

df=dt[dt$probeID%in%SigCpG, sampleSheet$Accession]

tdf=data.frame(t(df))
tdf$sample_ID=rownames(tdf)
tdf$Class=merge(tdf, sampleSheet, by.x="sample_ID", by.y="Accession")$Class
trainDT <- arrange(tdf, desc(Class))

write.table(trainDT, paste0(output_dir, "/Fig1_heatmap.tsv"), quote=F, row.names=F, sep="\t")