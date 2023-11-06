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
directory = "/path/to/data"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/fig/Fig1/data/B")

### Get DMPs
cpgList=read.table(paste0(directory, "/data/train/DMP.current/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)

DMP.all.GR=makeGRangesFromDataFrame(cpgList, keep.extra.columns=TRUE)
genome(DMP.all.GR)="hg38"

### Annotation
annots=c("hg38_genes_promoters", "hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_exons", "hg38_genes_introns", "hg38_genes_intergenic", "hg38_genes_intronexonboundaries", "hg38_genes_exonintronboundaries", "hg38_enhancers_fantom", "hg38_cpg_shores", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_inter", "hg38_cpgs")
annotations = build_annotations(genome = 'hg38', annotations = annots)
DMP.all.anno = annotate_regions(
  regions = DMP.all.GR,
  annotations = annotations,
  ignore.strand = FALSE,
  quiet = FALSE)

DMP.all.anno.df1 <- data.frame(DMP.all.anno)
dim(DMP.all.anno.df1)
table(DMP.all.anno.df1$annot.type)

DMP.all.anno.df <-  DMP.all.anno.df1[!duplicated(DMP.all.anno.df1[,c("probeID", "annot.type")]),]
dim(DMP.all.anno.df)
table(DMP.all.anno.df$annot.type)

### Split into relative CpG Island or Genes
## Relative CpG Island
DMP.cpg <- DMP.all.anno.df[grepl("hg38_cpg_", DMP.all.anno.df$annot.type),]
DMP.cpg$cpg.type <- "Non CpG islands"
DMP.cpg["cpg.type"][DMP.cpg["annot.type"] == "hg38_cpg_inter"] <- "Open sea"
DMP.cpg["cpg.type"][DMP.cpg["annot.type"] == "hg38_cpg_shores"] <- "CpG shores"
DMP.cpg["cpg.type"][DMP.cpg["annot.type"] == "hg38_cpg_islands"] <- "CpG islands"
DMP.cpg["cpg.type"][DMP.cpg["annot.type"] == "hg38_cpg_shelves"] <- "CpG shelves"
table(DMP.cpg$cpg.type)
fwrite(DMP.cpg, paste0(output_dir,"/annotation_cpgsites.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

## Genes
DMP.gene <- DMP.all.anno.df[!grepl("hg38_cpg_", DMP.all.anno.df$annot.type),]
DMP.gene$gene.type <- "Gene body"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_introns"] <- "Introns"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_intronexonboundaries"] <- "Intronic boundary"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_exonintronboundaries"] <- "Exonic boundary"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_intergenic"] <- "Intergenic"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_exons"] <- "Exons"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_enhancers_fantom"] <- "Enhancer (Fantom)"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_5UTRs"] <- "5UTR"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_promoters"] <- "Promoters"
DMP.gene["gene.type"][DMP.gene["annot.type"] == "hg38_genes_3UTRs"] <- "3UTR"
table(DMP.gene$gene.type)
fwrite(DMP.gene, paste0(output_dir,"/annotation_genes.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)