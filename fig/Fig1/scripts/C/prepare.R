library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
library(missMethyl)
library(data.table)
library(ggplot2)
set.seed(123)

### Setup directory
directory = "/path/to/data"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/fig/Fig1/data/C")

cpgList=read.table(paste0(directory, "/data/train/DMP.current/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)$probeID

ann=getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
allcpgs=rownames(ann)

# GO
GO_term=gometh(
  cpgList,
  all.cpg = NULL,
  collection = "GO",
  array.type = "EPIC",
  plot.bias = FALSE,
  prior.prob = TRUE,
  anno = ann,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  genomic.features = c("ALL"),
  sig.genes = TRUE
)

write.table(GO_term, paste0(output_dir,"/pathway_GO.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

