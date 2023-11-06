library(dplyr)
library(wesanderson)
library(ggplot2)
library(data.table)
library(annotatr)
library(GenomicRanges)
library(minfi)
library(stringr)
library(writexl)
library(RColorBrewer)

set.seed(123)

### Setup directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory = "/path/to/data_dir"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/data/train/DMP/TB-HC")

### Get DMPs
cpgList=read.table(paste0(directory, "/data/train/DMP/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)

DMP.all.GR=makeGRangesFromDataFrame(cpgList, keep.extra.columns=TRUE)
genome(DMP.all.GR)="hg38"

## Merge
# DMP.all.GR <- unlist(as(list(DMP.backgr.GR, DMP.hyper.GR, DMP.hypo.GR), "GRangesList"))

### Annotation        
# Genic
annots=c("hg38_genes_promoters", "hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_exons", "hg38_genes_introns", "hg38_genes_intergenic", "hg38_genes_intronexonboundaries", "hg38_genes_exonintronboundaries", "hg38_enhancers_fantom")

annotations=build_annotations(genome = 'hg38', annotations = annots)
DMP.all.anno=annotate_regions(
  regions=DMP.all.GR,
  annotations = annotations,
  ignore.strand = FALSE,
  quiet = FALSE)

DMP.all.anno.df1 <- data.frame(DMP.all.anno)

dm_annsum = as.data.frame(summarize_annotations(
    annotated_regions = DMP.all.anno.df1,
    quiet = TRUE))
dm_annsum$percent=round(100*(dm_annsum$n/sum(dm_annsum$n)),2)
dm_annsum$annot.type=c("Enhancer (Fantom)", "3UTR", "5UTR", "Exon boundary", "Exons", "Intergenics", "Intron boundary", "Introns", "Promoters")

color <- brewer.pal(nrow(dm_annsum), "Set3")
pdf(file = paste0(output_dir, "/GenicAnnotations.pdf"), width=9)
pie_labels <- paste0(dm_annsum$annot.type, " = ", dm_annsum$percent, "%")
pie(dm_annsum$percent, main="Genic Annotations", labels = pie_labels, col = color)
dev.off()

# CpG
annots=c("hg38_cpg_shores", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_inter", "hg38_cpgs")
annotations=build_annotations(genome = 'hg38', annotations = annots)

DMP.all.anno=annotate_regions(
  regions=DMP.all.GR,
  annotations = annotations,
  ignore.strand = FALSE,
  quiet = FALSE)

DMP.all.anno.df1 <- data.frame(DMP.all.anno)

# Randomize the input regions
dm_annsum = as.data.frame(summarize_annotations(
    annotated_regions = DMP.all.anno.df1,
    quiet = TRUE))
dm_annsum$percent=round(100*(dm_annsum$n/sum(dm_annsum$n)),2)
dm_annsum$annot.type=c("Open sea", "CpG island", "CpG shelves", "CpG shores")

color <- brewer.pal(nrow(dm_annsum), "Set3")
pdf(file = paste0(output_dir, "/CpG.Annotations.pdf"), width=9)
pie_labels <- paste0(dm_annsum$annot.type, " = ", dm_annsum$percent, "%")
pie(dm_annsum$percent, main="CpG Annotations", labels = pie_labels, col = color)
dev.off()