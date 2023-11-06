### data
awk '{OFS="\t"} NR>1&&$1!="NA" {print $1, $2, $3, $5}' /path/to/probe_annotation/hg38/EPIC_full.tsv | bedtools intersect -a stdin -b <(echo -e "chr7\t17296633\t17299500") | cut -f4 > /path/to/data/fig/Fig1/data/D/adjection_CpG.list

Rscript -e '
library(data.table)
directory = "/path/to/data"

input_dir=paste0(directory, "/fig/Fig1/data/D")
output_dir=paste0(directory, "/fig/Fig1/data/D")

TB=as.data.frame(fread(paste0(directory, "/data/dataInput/bigTable.mergedAll.tsv")))
CpG=data.frame(fread(paste0(input_dir, "/adjection_CpG.list"), header=F))
SampleSheet=read.table(paste0(directory, "/config/train.tsv"), sep="\t", header=T)

tmp=TB[,c("probeID", SampleSheet[,1])]

df=merge(CpG, tmp, by.x="V1", by.y="probeID", all.x=T)

df[is.na(df)] = 0
names(df)[1]="probeID"

write.table(df, paste0(output_dir, "Bval.tsv"), quote=F, row.names=F, sep="\t")
'