library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)

directory = "/path/to/data"

# TB dataset
TB=as.data.frame(fread(paste0(directory, "/data/dataInput", "/bigTable.mergedAll.tsv")))
# TB_Sheet1=read.table(paste0(directory, "/", "/config/train.tsv"),sep="\t", header=T)
# TB_Sheet2=read.table(paste0(directory, "/", "/config/test.tsv"),sep="\t", header=T)

# TB_Sheet=rbind(TB_Sheet1, TB_Sheet2)
# write.table(TB_Sheet, directory, "/fig/Fig4/data/BCG_Vaccination.tsv", sep="\t", quote=F, row.names=F)

pheno=read.table(directory, "/fig/Fig4/data/BCG_Vaccination.tsv", sep="\t", header=T)

set1 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/BMA.txt"), header=T)$probe_ID[-1]
set2 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/RF_min.txt"), header=T)$name
set3 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/PLS_min.txt"), header=T)$name
set4 <- read.table(paste0(directory, "/data/train/DMP/TB-HC/LASSO.txt"), header=T)$probeID

SigCpG=union(union(union(set1, set2), set3), set4)
supp=grep("cg21184174|cg10140638|cg10453758|cg19616230|cg00782174|cg15705999", SigCpG, invert=T, value=T)

TB_dt=TB[TB$probeID%in%supp, names(TB)%in%c("probeID", pheno$Accession)]

tmp=melt(TB_dt)
dat=merge(tmp, pheno, by.x="variable", by.y="Accession")

# ggplot(dat, aes(x=Class, y=value, fill=Class)) + 
#     geom_violin() +
#     facet_wrap(~probeID, scale="free")

# plot selected 
label3 <- list(c("TB", "BCG vaccination"), c("BCG vaccination", "HC") )

ggplot(dat[dat$probeID%in%supp,], aes(x=Class, y=value, fill=Class)) + 
geom_violin()  +
stat_compare_means(comparisons = label3, size=5) +
facet_wrap(~probeID, scale="free", ncol = 2) + ylab("Beta value") +
theme(text = element_text(size = 30), legend.position = "none") 
ggsave(directory, "/fig/FigS4/figure/HC_BCG.selected.pdf", width=15, height=20)
ggsave(directory, "/fig/FigS4/figure/HC_BCG.selected.png", width=15, height=20)
