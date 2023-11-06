# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install(c("minfi", "DMRcate", "missMethyl"))

library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(Gviz)
library(DMRcate)
library(stringr)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(ggfortify)
# library(config)

data_dir="/path/to/data_dir"

### load config
Pvalue=0.01
C=2
lambda=100
fdr=10^-8
# case=config$case
genome1="hg19"
genome2="hg38"
array="EPIC"
anno_dir="/path/to/probe_annotation"

### prepare design.csv
targets=read.table(paste0(data_dir, "/config/train.tsv"), header=T, sep="\t", check.names=F, stringsAsFactors=F)
# make sure no space in group's names
targets$base=str_replace_all(targets[,"Class"], " ", ".")

designTSV=targets[,c("Accession", "base")]
names(designTSV)=c("title", "groups")
write.table(designTSV, paste0(data_dir, "/config/design_train_data.tsv"), sep="\t", quote=F, row.names=F)

### DMPs - DMRs calling

bVals=read.table(paste0(data_dir, "/data/dataInput/bigTable.mergedAll.tsv"), header=T, sep=" ", check.names=F, stringsAsFactors=F)
# possible comparisons
pair=combn(unique(designTSV$groups),2)

gr=pair[,1]
print(gr)

g1=gr[1]
g2=gr[2]

SubGr=designTSV[designTSV$groups%in%gr,]
sN=SubGr$title
rownames(bVals)=bVals$probeID
dt=bVals[,sN]

# matrix
type=ifelse(grepl(g1, SubGr$groups), "case", "control")
type=relevel(factor(type), "control")
design=model.matrix(~type)
colnames(design)=gsub("type", "", colnames(design))

# dt=dt[complete.cases(dt),]
df=dt[complete.cases(dt),]
df=df/100

myAnnotation=cpg.annotate("array", as.matrix(df), arraytype=array, analysis.type="differential", design=design, coef="case", fdr=fdr, what="Beta")
cpSites=as.data.frame(myAnnotation@ranges)

# fit <- lmFit(as.matrix(df), design=design)
# fit <- eBayes(fit)
# topTable(fit)

DMP.query=cpSites[cpSites$is.sig==TRUE,]

DMR=paste0(data_dir, "/data/train/DMR")
if (!dir.exists(DMR)) {dir.create(DMR)}
path=paste0(DMR, "/", paste0(g1, "-", g2))
if (!dir.exists(path)) {dir.create(path, recursive=T)}

dmrcate.res=dmrcate(myAnnotation, lambda=lambda, C=C, min.cpgs=2)
DMRquery=length(dmrcate.res@coord)

saveRDS(dmrcate.res, file = paste0(path, "/DMR.rds"))
### extract ranges
# genome1
dmrcate.ranges=extractRanges(dmrcate.res, genome=genome1)
write.table(dmrcate.ranges, paste0(path, "/DMR.", genome1, ".tsv"), quote=F, row.names=F, sep="\t")

# # genome2
dmrcate.ranges=extractRanges(dmrcate.res, genome=genome2)
write.table(data.frame(dmrcate.ranges), paste0(path, "/DMR.", genome2, ".tsv"), quote=F, row.names=F, sep="\t")

### DMP

print("DMPs calling")
DMPs=cpSites[cpSites$is.sig==TRUE,]
DMPs$probeID=rownames(DMPs)

DMP=paste0(data_dir, "/data/train/DMP")
if (!dir.exists(DMP)) {dir.create(DMP, recursive=T)}

path=paste0(data_dir, "/data/train/DMP/", paste0(g1, "-", g2))
if (!dir.exists(path)) {dir.create(path, recursive=T)}


DMPs$logFC=merge(DMPs, FoldChange, by=0)$logFC

write.table(DMPs, paste0(path, "/DMPs.", "hg38.tsv"), quote=F, sep="\t", row.names=F)