library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)

# TB dataset
TB=as.data.frame(fread(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/data/dataInput", "/bigTable.mergedAll.tsv")))
# TB_Sheet1=read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/", "/config/train.tsv"),sep="\t", header=T)
# TB_Sheet2=read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/", "/config/test.tsv"),sep="\t", header=T)

# TB_Sheet=rbind(TB_Sheet1, TB_Sheet2)[,c(1,2,3)]
# TB_Sheet=TB_Sheet[TB_Sheet$Class=="TB",]
TB_Sheet=read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/", "/config/totalTB.tsv"), sep="\t", header=T)[,c(1,7)]

SigCpG=c("cg00782174", "cg01737507",  "cg05501357",  "cg06270401",  "cg06497752",  "cg09418321", "cg10140638",  "cg10453758", "cg11166252",  "cg14095850",  "cg15705999",  "cg18644543",  "cg19616230",  "cg20098659",  "cg21184174",  "cg22082462",  "cg22381196", "cg23181133")

TB_dt=TB[TB$probeID%in%SigCpG, names(TB)%in%c("probeID", TB_Sheet$Accession)]

# HIV dataset
HIV=as.data.frame(fread(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/data/GSE53840/methProbes", "/bigTable.hg38.tsv")))[,-c(1:4)]

HIV_Sheet=read.table(paste0("/media/thong/sda/RCID/CoHoai/Tuberculosis/data/GSE53840/raw", "/SampleSheet.tsv"),sep="\t", header=T)

df=merge(TB_dt, HIV, by="probeID")
sheet=rbind(TB_Sheet, HIV_Sheet)

ldf=melt(df)
dat=merge(ldf, sheet, by.x="variable", by.y="Accession")

ggplot(dat, aes(x=Class, y=value, fill=Class)) + 
    geom_violin() +
    facet_wrap(~probeID, scale="free")

label3 <- list(c("HC", "TB"), c("TB", "HIV"), c("TB", "TB+HIV"), c("HC", "TB+HIV") )

ggplot(dat, aes(x=Class, y=value, fill=Class)) + 
geom_violin()  +
stat_compare_means(comparisons = label3) +
facet_wrap(~probeID, scale="free")
ggsave("/media/thong/sda/RCID/CoHoai/Tuberculosis/TB_HIV.plot.pdf", =25)
