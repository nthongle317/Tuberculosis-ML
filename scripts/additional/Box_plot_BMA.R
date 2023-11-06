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
set.seed(123)

directory = "/path/to/data_dir"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/data/test/Model/BMA")

dt=as.data.frame(fread(paste0(input_dir, "/bigTable.mergedAll.tsv")))
row.names(dt)=dt$probeID
SigCpG=read.table(paste0(directory, "/data/train/DMP/TB-HC/BMA.txt"), sep="\t", header=T)$probe_ID

# validation data
ValSample=read.table(paste0(directory, "/config/test.tsv"),sep="\t", header=T)
df=dt[dt$probeID%in%SigCpG, ValSample$Accession]
tdf=data.frame(t(df))
tdf$sample_ID=rownames(tdf)
tdf$Class=merge(tdf, ValSample, by.x="sample_ID", by.y="Accession")$Class
valDT=arrange(tdf, desc(Class))
valDT=valDT[,-grep("sample_ID", names(valDT))]

test_plot=melt(valDT)

scaleFUN <- function(x) sprintf("%.1f", x)

### Plot

pdf(file=file.path(paste0(output_dir, "/", "boxplot_2classs_BMA.pdf")), width = 8, height = 7)

ggplot(test_plot, aes(x = Class, y = value, fill = Class)) +
geom_boxplot(width = 0.15, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable) +
stat_compare_means(label.y = 105,size = 9) +
scale_y_continuous(labels=scaleFUN, limits = c(0, 120)) +  
ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") #+
# stat_n_text(size = 9, y.pos = 105)

dev.off()

# violin only
ggplot(test_plot, aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.15, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable) +
stat_compare_means(label.y = 105,size = 9) +
scale_y_continuous(labels=scaleFUN, limits = c(0, 120)) +  
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") #+
# stat_n_text(size = 9, y.pos = 105)

# violin only
ggplot(test_plot, aes(x = Class, y = value, fill = Class)) +
geom_violin(width = 0.15, color="black",lwd = 1, outlier.color = NA) + facet_wrap(~variable) +
stat_compare_means(label.y = 105,size = 9) +
scale_y_continuous(labels=scaleFUN, limits = c(0, 120)) +  
# ggdist::stat_halfeye(aes(fill = Class), adjust = .5, width = .3, justification = -.6,.width = 0, point_colour = NA) + 
# gghalves::geom_half_point(aes(color = Class),side = "l", range_scale = .3, size = 3, alpha = .6) +
labs(fill = "",color = "") + xlab("") + ylab(paste0("Beta value")) +
scale_fill_manual(values=c("#CC79A7","#FC4E07")) +
scale_color_manual(values = c("#CC79A7", "#FC4E07")) +
theme(text = element_text(size = 30), legend.position = "none") #+
