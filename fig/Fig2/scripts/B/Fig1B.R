library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(RColorBrewer)

set.seed(123)

### Setup directory
directory = "/path/to/data"

input_dir=paste0(directory, "/fig/Fig1/data/B")
output_dir=paste0(directory, "/fig/Fig1/figure/B")

### Get CpG sites and genes annotation
## CpG sites annotation
DMP.cpg <- as.data.frame(fread(paste0(input_dir, "/annotation_cpgsites.tsv")))

df=data.frame(table(DMP.cpg$cpg.type))
df$percent=round(100*(df$Freq/sum(df$Freq)),2)

# plot
color <- brewer.pal(nrow(df), "Set3")
pdf(file = paste0(output_dir, "/CpGs.Annotations.pdf"), width=9)
pie_labels <- paste0(df$Var1, " = ", df$percent, "%")
pie(df$percent, main="CpG Annotations", labels = pie_labels, col = color)
dev.off()

png(file = paste0(output_dir, "/CpGs.Annotations.png"), width=700, height=600)
pie_labels <- paste0(df$Var1, " = ", df$percent, "%")
pie(df$percent, main="CpG Annotations", labels = pie_labels, col = color)
dev.off()


## Genic annotation
DMP.gene <- as.data.frame(fread(paste0(input_dir, "/annotation_genes.tsv")))

df=data.frame(table(DMP.gene$gene.type))
df$percent=round(100*(df$Freq/sum(df$Freq)),2)

# plot
color <- brewer.pal(nrow(df), "Set3")
pdf(file = paste0(output_dir, "/Genic.Annotations.pdf"), width=9)
pie_labels <- paste0(df$Var1, " = ", df$percent, "%")
pie(df$percent, main="Genic Annotations", labels = pie_labels, col = color)
dev.off()

png(file = paste0(output_dir, "/Genic.Annotations.png"), width=700, height=600)
pie_labels <- paste0(df$Var1, " = ", df$percent, "%")
pie(df$percent, main="Genic Annotations", labels = pie_labels, col = color)
dev.off()