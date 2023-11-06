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

input_dir=paste0(directory, "/fig/FigS2/data/S2")
output_dir=paste0(directory, "/fig/FigS2/figure/S2")

dt=data.frame(fread(paste0(input_dir, "/boxplotData.tsv")))
FDR=data.frame(fread(paste0(input_dir, "/Limma_indFDR.tsv")))
FDR$ind.fdr=signif(FDR$ind.fdr, digits = 3)

forPlot=grep("cg23181133|cg20098659|cg10453758|cg19616230|cg21184174|cg10140638", FDR$probeID, value=T, invert=T)
df=dt[dt$variable%in%forPlot,]


# plot all
scaleFUN <- function(x) sprintf("%.1f", x)
ggplot(df, aes(x = Class, y = value, fill = Class)) +
geom_violin(width=0.3, color="black",lwd = 1) + facet_wrap(~variable, scale='free', ncol=4) +
stat_compare_means(comparisons = list(c("HC", "TB")), size=5) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7|#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") 
ggsave(file.path(paste0(output_dir, "/", "figureS2.pdf")), width=10, height=8)
