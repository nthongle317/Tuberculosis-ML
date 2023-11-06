library(ggplot2)
library(data.table)
library(annotatr)
library(GenomicRanges)
library(minfi)
library(stringr)
library(writexl)
library(RColorBrewer)
library(dplyr)

set.seed(123)

### Setup directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory = "/path/to/data"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/fig/Fig1/data/D")

### Get data
cpgList=read.table(paste0(directory, "/data/train/DMP.current/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)
cpgList$start=cpgList$start-500
cpgList$end=cpgList$end+500
cpgList$width=cpgList$end-cpgList$start

DMP.all.GR=makeGRangesFromDataFrame(cpgList, keep.extra.columns=TRUE)
DMP.all.GR$dmrID <- paste0("DMR_", sprintf("%06.0f", 1:length(DMP.all.GR)))

genome(DMP.all.GR)="hg38"

### Annotation
annots=c("hg38_genes_promoters", "hg38_genes_3UTRs", "hg38_genes_5UTRs", "hg38_genes_exons", "hg38_genes_introns", "hg38_genes_intergenic", "hg38_genes_intronexonboundaries", "hg38_genes_exonintronboundaries", "hg38_enhancers_fantom", "hg38_cpg_shores", "hg38_cpg_islands", "hg38_cpg_shelves", "hg38_cpg_inter", "hg38_cpgs")
annotations = build_annotations(genome = 'hg38', annotations = annots)
DMP.all.anno = annotate_regions(
  regions = DMP.all.GR,
  annotations = annotations,
  ignore.strand = FALSE,
  quiet = FALSE)

DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_cpg_inter"] <- "Open sea"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_cpg_shores"] <- "CpG shores"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_cpg_islands"] <- "CpG islands"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_cpg_shelves"] <- "CpG shelves"
DMP.all.anno$annot$type[DMP.all.anno$annot$type== "hg38_genes_introns"] <- "Introns"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_genes_intronexonboundaries"] <- "Intronic boundary"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_genes_exonintronboundaries"] <- "Exonic boundary"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_genes_intergenic"] <- "Intergenic"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_genes_exons"] <- "Exons"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_enhancers_fantom"] <- "Enhancer (Fantom)"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_genes_5UTRs"] <- "5UTR"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_genes_promoters"] <- "Promoters"
DMP.all.anno$annot$type[DMP.all.anno$annot$type=="hg38_genes_3UTRs"] <- "3UTR"

## Filter 2: 
# CpG islands DMRs only
dmr.island <- DMP.all.anno[DMP.all.anno$annot$type=="CpG islands", "dmrID"]
# Promoter DMRs only
dmr.promoter <- DMP.all.anno[DMP.all.anno$annot$type=="Promoters", "dmrID"]
# select markers
# markers=c("cg00782174", "cg01737507", "cg05501357", "cg06270401", "cg06497752", "cg09418321", "cg10140638", "cg10453758", "cg11166252", "cg14095850", "cg15705999", "cg18644543", "cg19616230", "cg20098659", "cg21184174", "cg22082462", "cg22381196", "cg23181133")
# dmr_markers <- DMP.all.anno[DMP.all.anno$probeID%in%markers & , "dmrID"]
# select Enhancer
# dmr_markers <- DMP.all.anno[DMP.all.anno$probeID%in%markers, "dmrID"]

# Intersection
dmr.interested <- intersect(dmr.promoter$dmrID, dmr_markers$dmrID)
sig.dmr <- DMP.all.anno[DMP.all.anno$dmrID %in% dmr.interested]
length(unique(sig.dmr@ranges))
saveRDS(sig.dmr, file = paste0(output_dir, "/DMR_forplot.rds"))

### Choose my ranges
sig.start <- ifelse(sig.dmr@ranges@start < sig.dmr$annot@ranges@start,
                    sig.dmr@ranges@start,
                    sig.dmr$annot@ranges@start)
sig.end <- ifelse((sig.dmr@ranges@start+sig.dmr@ranges@width) >
                    (sig.dmr$annot@ranges@start+sig.dmr$annot@ranges@width),
                  (sig.dmr@ranges@start+sig.dmr@ranges@width),
                  (sig.dmr$annot@ranges@start+sig.dmr$annot@ranges@width))
sig.ranges <- GRanges(sig.dmr@seqnames, IRanges(sig.start, sig.end))

### Get phenotype
pheno_type <- as.data.frame(fread(paste0(directory, "/config/train.tsv")))
pheno_type <- arrange(pheno_type, (Class))

### Get beta value
data=as.data.frame(fread(paste0(directory, "/data/dataInput/bigTable.mergedAll.tsv")))
rownames(data) <- data$probeID
data <- data[, c("probeID", pheno_type$Accession)]
dt <- data[,-1]

## Turn beta value matrix to Granges format
GRanges_DMP <- makeGenomicRatioSetFromMatrix(as.matrix(dt), array = "IlluminaHumanMethylationEPIC",
                                             annotation = "ilm10b5.hg38", mergeManifest = TRUE, what = "Beta")
CpGs <- getBeta(GRanges_DMP)
RSanno <- getAnnotation(GRanges_DMP)
RSanno <- RSanno[order(RSanno$chr, RSanno$pos), ]
CpGs <- CpGs[rownames(RSanno), ]
cpgs.ranges <- GRanges(RSanno$chr, IRanges(RSanno$pos, 
                                           RSanno$pos))
values(cpgs.ranges) <- CpGs

## Intersect probes and forplot DMRs
cpgs.ranges <- subsetByOverlaps(cpgs.ranges, sig.ranges)
saveRDS(cpgs.ranges, file = paste0(output_dir, "/DMP_forplot.rds"))
