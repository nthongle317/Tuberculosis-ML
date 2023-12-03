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

input_dir=paste0(directory, "/fig/Fig2/data/B")
output_dir=paste0(directory, "/fig/Fig2/figure/B")

dt=data.frame(fread(paste0(input_dir, "/boxplotData.tsv")))
FDR=data.frame(fread(paste0(input_dir, "/Limma_indFDR.tsv")))
FDR$ind.fdr=signif(FDR$ind.fdr, digits = 3)
scaleFUN <- function(x) sprintf("%.1f", x)

# plot all
for (cpg in FDR$probeID) {

ggplot(dt[dt$variable==cpg,], aes(x = Class, y = value, fill = Class)) +
geom_violin(width=0.5, color="black",lwd = 1) + facet_wrap(~variable, scale='free') +
stat_compare_means(comparisons = list(c("HC", "TB")), size=7) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") 
ggsave(file.path(paste0(output_dir, "/", cpg,  "_boxplot.pdf")), width=10)

}

# plot selected
cpg="cg23181133"
p1=ggplot(dt[dt$variable==cpg,], aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.5, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable, scale='free') +
stat_compare_means(comparisons = list(c("HC", "TB")), size=7) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") 

cpg="cg20098659"
p2=ggplot(dt[dt$variable==cpg,], aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.5, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable, scale='free') +
stat_compare_means(comparisons = list(c("HC", "TB")), size=7) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") 

cpg="cg10453758"
p3=ggplot(dt[dt$variable==cpg,], aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.5, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable, scale='free') +
stat_compare_means(comparisons = list(c("HC", "TB")), size=7) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") 

cpg="cg19616230"
p4=ggplot(dt[dt$variable==cpg,], aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.5, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable, scale='free') +
stat_compare_means(comparisons = list(c("HC", "TB")), size=7) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none")

cpg="cg21184174"
p5=ggplot(dt[dt$variable==cpg,], aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.5, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable, scale='free') +
stat_compare_means(comparisons = list(c("HC", "TB")), size=7) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none")

cpg="cg10140638"
p6=ggplot(dt[dt$variable==cpg,], aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.5, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable, scale='free') +
stat_compare_means(comparisons = list(c("HC", "TB")), size=7) +
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none")


ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)
ggsave(file.path(paste0(output_dir, "/", "boxplot_selected.pdf")), width=12, height=20)
ggsave(file.path(paste0(output_dir, "/", "boxplot_selected.png")), width=12, height=20)