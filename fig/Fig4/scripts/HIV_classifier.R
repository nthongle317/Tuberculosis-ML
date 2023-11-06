library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)

directory = "/path/to/data"

# TB dataset
TB=as.data.frame(fread(paste0(directory, "/data/dataInput", "/bigTable.mergedAll.tsv")))
# TB_Sheet1=read.table(paste0(directory, "/", "/config/train.tsv"),sep="\t", header=T)
# TB_Sheet2=read.table(paste0(directory, "/", "/config/test.tsv"),sep="\t", header=T)

# TB_Sheet=rbind(TB_Sheet1, TB_Sheet2)[,c(1,2,3)]
# TB_Sheet=TB_Sheet[TB_Sheet$Class=="TB",]
TB_Sheet=read.table(paste0(directory, "/", "/config/totalTB.tsv"), sep="\t", header=T)[,c(1,7)]

set1 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/BMA.txt"), header=T)$probe_ID[-1]
set2 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/RF_min.txt"), header=T)$name
set3 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/PLS_min.txt"), header=T)$name
set4 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/LASSO.txt"), header=T)$probeID

SigCpG=union(union(union(set1, set2), set3), set4)

TB_dt=TB[TB$probeID%in%SigCpG, names(TB)%in%c("probeID", TB_Sheet$Accession)]

# HIV dataset
HIV=as.data.frame(fread(paste0(directory, "/data/GSE53840/methProbes", "/bigTable.hg38.tsv")))[,-c(1:4)]

HIV_Sheet=read.table(paste0(directory, "/data/GSE53840/raw", "/SampleSheet.tsv"),sep="\t", header=T)

df=merge(TB_dt, HIV, by="probeID")
sheet=rbind(TB_Sheet, HIV_Sheet)

ldf=melt(df)
dat=merge(ldf, sheet, by.x="variable", by.y="Accession")

ggplot(dat, aes(x=Class, y=value, fill=Class)) + 
    geom_violin() +
    facet_wrap(~probeID, scale="free")

label3 <- list(c("TB", "HIV"), c("TB", "TB+HIV"), c("HC", "TB+HIV") )

ggplot(dat, aes(x=Class, y=value, fill=Class)) + 
geom_violin()  +
stat_compare_means(comparisons = label3, size=7) +
facet_wrap(~probeID, scale="free") +
theme(text = element_text(size = 30), legend.position = "none") 
ggsave(directory, "/fig/Fig4/figure/TB_HIV.plot.pdf", width=25, height=25)


# plot selected 
label3 <- list(c("TB", "HIV"), c("TB", "TB+HIV"), c("HC", "TB+HIV") )

ggplot(dat[dat$probeID%in%c("cg00782174", "cg10140638", "cg10453758", "cg19616230", "cg20098659", "cg23181133"),], aes(x=Class, y=value, fill=Class)) + 
geom_violin()  +
stat_compare_means(comparisons = label3, size=7) +
facet_wrap(~probeID, scale="free", ncol = 2) + ylab("Beta value") +
theme(text = element_text(size = 30), legend.position = "none") 
ggsave(directory, "/fig/Fig4/figure/TB_HIV.selected.pdf", width=15, height=20)
ggsave(directory, "/fig/Fig4/figure/TB_HIV.selected.png", width=15, height=20)
