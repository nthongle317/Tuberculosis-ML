library(data.table)
library(dplyr)
library(caret)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

set.seed(123)

### Setup directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory = "/path/to/data_dir"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/data/train/DMP/TB-HC")

### get beta matrix
dt <- as.data.frame(fread(paste0(input_dir, "/bigTable.mergedAll.tsv")))
row.names(dt)=dt$probeID
sampleSheet=read.table(paste0(directory, "/config/train.tsv"),sep="\t", header=T)
cpgList=read.table(paste0(directory, "/data/train/DMP/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)

SigCpG=cpgList[abs(cpgList$diff)>=0.1,"probeID"]

cpgList=cpgList[cpgList$probeID%in%SigCpG, ]

df=dt[dt$probeID%in%SigCpG, sampleSheet$Accession]

tdf=data.frame(t(df))
tdf$sample_ID=rownames(tdf)
tdf$Class=merge(tdf, sampleSheet, by.x="sample_ID", by.y="Accession")$Class
trainDT <- arrange(tdf, desc(Class))

t.train <- t(trainDT[,SigCpG])

### Set color and annotation
# gradient color for beta value
col_fun=colorRamp2(c(0, 50, 100), c("darkblue", "white", "darkred"))
# col_fun(seq(-3, 3))

# Split column
no.TB=NROW(trainDT[trainDT$Class == "TB",])
# no.Mild=NROW(trainDT[trainDT$status == "Mild",])
no.HC=NROW(trainDT[trainDT$Class == "HC",])
split.class=c(rep("TB", each=no.TB),
              # rep("Mild", each=no.Mild),
              rep("HC", each=no.HC))

# Split row
no.hyper=NROW(cpgList[cpgList$diff <= 0,])
no.hypo=NROW(cpgList[cpgList$diff >= 0,])
split.trend=c(rep("Hypo DMPs", each=no.hypo),
              rep("Hyper DMPs", each=no.hyper))
cpgList$trend=ifelse(cpgList$diff <= 0, "Hyper DMPs", "Hypo DMPs")
cpgList <- arrange(cpgList, desc(trend))

# class label
classcolor=list(bar=c("TB"="#FC4E07", "HC"="#00AFBB"))

ha.col=HeatmapAnnotation(bar=trainDT[,"Class"],
                         show_annotation_name=FALSE,
                         show_legend=FALSE,
                         border = FALSE,
                         col=classcolor)
draw(ha.col)

# trend label
row.anno.df=data.frame(Trend=cpgList[,"trend"])
ha.row = rowAnnotation(df = row.anno.df,
                       col = list(Trend = c("Hypo DMPs"="#D8E9A8", "Hyper DMPs"="purple")),
                       show_annotation_name=FALSE,
                       show_legend=TRUE,
                       border = FALSE,
                       width = unit(1, "cm"))
draw(ha.row)


### Plot heatmap
heatmap_plot <- Heatmap(t.train, name="Methylation",
                        col=col_fun,
                        top_annotation=ha.col,
                        left_annotation = ha.row, 
                        column_split=split.class,
                        row_split=split.trend,
                        show_column_names=FALSE, 
                        show_row_names=FALSE,
                        cluster_columns=FALSE,
                        cluster_rows=TRUE,
                        border = FALSE)


png(file=file.path(paste0(output_dir, "/heatmap_0.1.png")), width=800, height=600)
print(heatmap_plot)
dev.off()

pdf(file=file.path(paste0(output_dir, "/heatmap_0.1.pdf")), width=8, height=7)
print(heatmap_plot)
dev.off()
# }

