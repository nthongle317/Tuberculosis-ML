library(dplyr)
library(wesanderson)
library(ggplot2)
library(data.table)
library(annotatr)
library(GenomicRanges)
library(Gviz)
library(minfi)
library(stringr)

set.seed(123)
genome = "hg38"

### Setup directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory = "/path/to/data"

input_dir=paste0(directory, "/fig/Fig1/data/D")
output_dir=paste0(directory, "/fig/Fig1/figure/D")

### Get phenotype
pheno_type <- as.data.frame(fread(paste0(directory, "/config/train.tsv")))
pheno_type <- arrange(pheno_type, (Class))

### Get forplot DMPs data
DMPs_anno <- readRDS(paste0(input_dir, "/DMP_forplot.rds"))

### Get forplot DMRs data
DMRs_anno <- readRDS(paste0(input_dir, "/DMR_forplot.rds"))
DMRs_anno$annot$cpg <- ifelse(grepl("CpG", DMRs_anno$annot$type),
                              DMRs_anno$annot$type, NA)

DMRs_anno$annot$genes <- ifelse(!grepl("CpG", DMRs_anno$annot$type),
                                DMRs_anno$annot$type, NA)

DMRs <- DMRs_anno[!duplicated(DMRs_anno$dmrID)]


### Plot by forloop
DMR_ID <- unique(DMRs_anno$dmrID)

for (i in seq(DMR_ID)) {
  i_dmr <- DMR_ID[i]
  print(i_dmr)
  
  ## Prepare something
  dmr_plot <- DMRs_anno[DMRs_anno$dmrID == i_dmr]
  chr = as.character(unique(dmr_plot@seqnames))
  start=unique(dmr_plot@ranges@start)
  width=unique(dmr_plot@ranges@width)
  end = unique(dmr_plot@ranges@start + dmr_plot@ranges@width)
  print(paste0(i_dmr, " located in ", chr, " (From ", start, " to ", end, ")"))
  
  ## Full tracks plot
  # Chromosome track
  itrack <- IdeogramTrack(genome = genome, chromosome = chr, fontcolor = "black", cex=2)
  
  # Genomic track
  gatrack <- GenomeAxisTrack(cex=1.5)
  
  # DMR track
  DMRtrack <- AnnotationTrack(dmr_plot[1], name = "DMR",
                              showFeatureId = FALSE,
                              fill = "purple",
                              col = "black", 
                              background.title = "#FFFEDB",
                              col.title = "black",
                              rotation.title=0,
                              cex.title=1.5,
                              background.panel = "#FFFEDB",
                              col.bosrder.title="#FFFEDB")
  plotTracks(DMRtrack)
  
  # CpG Island track
  dmr_cpg <- dmr_plot$annot[!is.na(dmr_plot$annot$cpg)]
  island_names <- unique(dmr_cpg$cpg)
  print(paste0(i_dmr, " overlap ", island_names))
  
  islandtrack <- AnnotationTrack(dmr_cpg[1],
                                 name = "CpG island",
                                 fill = "darkgreen",
                                 showFeatureId = FALSE,
                                 col = "black", 
                                 background.title = "#FFFEDB",
                                 col.title = "black",
                                 rotation.title=0,
                                 cex.title=1.5,
                                 background.panel = "#FFFEDB",
                                 col.border.title="#FFFEDB")
  
  plotTracks(islandtrack)
  
  # Gene track
  dmr_gene <- dmr_plot$annot[!is.na(dmr_plot$annot$genes)]
  dmr_gene <- dmr_plot$annot[!is.na(dmr_plot$annot$symbol)]
  dmr_gene <- dmr_gene[!duplicated(dmr_gene$genes),]
  print(paste0(i_dmr, " overlap gene element "))
  print(table(dmr_gene$genes, dmr_gene$symbol))
  genes_names <- paste0(unique(dmr_gene$symbol), " gene")
  dmr_gene2 <- as.data.frame(dmr_gene)
  
  # refGenes <- UcscTrack(genome = "hg38", chromosome = unique(dmr_gene2$seqnames),
  #                     track = "xenoRefGene", from = min(dmr_gene2$start), to = max(dmr_gene2$end),
  #                     trackType = "GeneRegionTrack", 
  #                     rstarts = "exonStarts", rends = "exonEnds", 
  #                     gene = "name",  symbol = "name2", 
  #                     transcript = "name", strand = "strand",
  #                     fill = "#8282d2", stacking = "dense", 
  #                     name = genes_names)


  if (length(genes_names) == 1) {
    dmr_gene2$symbol <- dmr_gene2$id
    genestrack <- GeneRegionTrack(dmr_gene2,
                                  fill = "salmon",
                                  chromosome=chr,
                                  genome=genome,
                                  name=genes_names,
                                  showId = TRUE,
                                  groupAnnotation="id",
                                  col = "black", 
                                  background.title = "transparent",
                                  col.title = "black",
                                  rotation.title=0,
                                  cex.title=1.5,
                                  just.group="above",
                                  cex.group=1)
  }
  
  
  if (length(genes_names) > 1) {
    dmr_gene2$symbol <- paste0(dmr_gene2$symbol, " (", dmr_gene2$id, ")")
    genestrack <- GeneRegionTrack(dmr_gene2,
                                  fill = "salmon",
                                  chromosome=chr,
                                  genome=genome,
                                  name="Genes",
                                  showId = TRUE,
                                  groupAnnotation="id",
                                  col = "black", 
                                  background.title = "transparent",
                                  col.title = "black",
                                  rotation.title=0,
                                  cex.title=1.5,
                                  just.group="above",
                                  cex.group=1)
  }
  
  
  plotTracks(genestrack)
  
  # DMPs track
  i_cpgs.ranges <- subsetByOverlaps(DMPs_anno, dmr_gene, maxgap=1000)
  # i_cpgs.ranges <- i_cpgs.ranges[-length(i_cpgs.ranges)]
  i_cpgs.ranges.df <- as.data.frame(i_cpgs.ranges)
  DMPtrack <- AnnotationTrack(i_cpgs.ranges,
                              chromosome = chr,
                              name = "DMPs",
                              fill = "green",
                              stacking = "dense",
                              showFeatureId = FALSE,
                              col = "transparent", 
                              background.title = "#FFFEDB",
                              col.title = "black",
                              rotation.title=0,
                              cex.title=1.5,
                              background.panel = "#FFFEDB",
                              col.border.title="#FFFEDB")
  plotTracks(DMPtrack)
  
  # Data track
  Nsample <- as.data.frame(table(pheno_type$Class))
  dtrack <- DataTrack(i_cpgs.ranges, 
                      chromosome = chr, 
                      genome = genome,
                      name = "DNA methylation", 
                      background.title = "#055350",
                      groups = c(rep("TB", Nsample[Nsample$Var1 == "TB", "Freq"]),
                                 rep("HC", Nsample[Nsample$Var1 == "HC", "Freq"])),
                      col=c("#00AFBB", "#F0E442"),
                      type = c("g", "a", "confint"),
                      legend = TRUE,
                      box.legend=TRUE,
                      cex.legend = 1.5,
                      cex = 1.5,
                      col.title = "white",
                      rotation.title=90,
                      cex.title=1.5,
                      cex.axis=1.5, lwd=3,
                      jitter.x=TRUE , factor=1, pch=3)
  plotTracks(dtrack)
  
  # Highlight track
  ht <- HighlightTrack(trackList = list(DMRtrack, DMPtrack, dtrack),
                       start = start-10, width = width+20,
                       chromosome = chr, col="black",
                       fill="transparent", lwd=0.5,
                       inBackground=FALSE)
  
  # Save it
  pdf(paste0(output_dir, "/", i_dmr, "_FullTrack.pdf"), width = 8, height = 7)
  plotTracks(list(itrack, gatrack, genestrack, islandtrack, ht))
  dev.off()
  
  # ## Zoom-in tracks plot
  # # Zoom-in DMPs track
  # i_cpgs.ranges2 <- subsetByOverlaps(DMPs_anno, dmr_plot[1])
  # print(paste0(i_dmr, " include ", length(i_cpgs.ranges2), " DMPs "))
  # i_cpgs.ranges2.df <- as.data.frame(i_cpgs.ranges2)
  # DMPtrack2 <- AnnotationTrack(i_cpgs.ranges2,
  #                              chromosome = chr,
  #                              name = "DMPs",
  #                              fill = "green",
  #                              stacking = "dense",
  #                              showFeatureId = FALSE,
  #                              col = "transparent", 
  #                              background.title = "#FFFEDB",
  #                              col.title = "black",
  #                              rotation.title=0,
  #                              cex.title=1.5,
  #                              background.panel = "#FFFEDB",
  #                              col.border.title="#FFFEDB")
  # plotTracks(DMPtrack2)
  
  # # Zoom-in data track
  # dtrack2 <- DataTrack(i_cpgs.ranges2, 
  #                      chromosome = chr, 
  #                      genome = genome,
  #                      name = "DNA methylation", 
  #                      background.title = "#055350",
  #                     groups = c(rep("HC", Nsample[Nsample$Var1 == "HC", "Freq"]),
  #                                rep("TB", Nsample[Nsample$Var1 == "TB", "Freq"])),
  #                      col=c("#00AFBB", "#F0E442"),
  #                      type = c("boxplot", "g", "a"), 
  #                      legend = TRUE, 
  #                      cex.legend = 1.5,
  #                      cex = 1.5,
  #                      col.title = "white",
  #                      rotation.title=90,
  #                      cex.title=1.5,
  #                      cex.axis=1.5, 
  #                      lwd=3,
  #                      box.legend=TRUE,
  #                      jitter.x=TRUE, 
  #                      factor=1, 
  #                      pch=3)
  # plotTracks(dtrack2)
  
  # # Save it
  # pdf(paste0(output_dir, "/", i_dmr, "_ZoomInTrack.pdf"), width = 8, height = 7)
  # plotTracks(list(itrack, gatrack, DMRtrack, DMPtrack2, dtrack2), 
  #            from = start-10, to = end+10)
  # dev.off()
  
}

