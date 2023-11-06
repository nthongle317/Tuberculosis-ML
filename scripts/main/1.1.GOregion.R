# library(devtools)
# install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest") 
# install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")

library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
library(missMethyl)
library(data.table)
library(ggplot2)
set.seed(123)

### Setup directory
directory = "/path/to/data_dir"

input_dir=paste0(directory, "/data/dataInput")
output_dir=paste0(directory, "/data/train/DMP/TB-HC")

cpgList=read.table(paste0(directory, "/data/train/DMP/TB-HC/DMPs.hg38.tsv"), sep="\t", header=T)$probeID

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

df=GO_term[GO_term$P.DE<=0.05,]

write.table(GO_term, paste0(output_dir,"/pathway_GO.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# plot
# Filtering
TFTerm_GO=df[grepl("tuberculosis", df$TERM, ignore.case=T) |
                               grepl("interleukin", df$TERM, ignore.case=T) |
                               grepl("chemokine", df$TERM, ignore.case=T) |
                               grepl("immune", df$TERM, ignore.case=T) | 
                               grepl("kappa", df$TERM, ignore.case=T) |
                               grepl("TNF", df$TERM, ignore.case=T) |
                               grepl("inflamatory", df$TERM, ignore.case=T),]
ggplot(TFTerm_GO,
       aes(x=reorder(TERM, DE, decreasing = FALSE), y=DE)) +
  geom_bar(aes(fill = P.DE), stat="identity") +
  scale_fill_gradientn(colours = rainbow(5)) +
  labs(x = '', y = "#of genes that are differentially methylated in the term", fill = "P.DE") +
  coord_flip() +
  facet_wrap(~ONTOLOGY) +
  theme(text = element_text(size = 20),
        axis.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold', size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(output_dir,"/region_pathway_GO.pdf"), width = 20, height = 10)


# KEGG
KEGG_term=gometh(
  cpgList,
  all.cpg = NULL,
  collection = "KEGG",
  array.type = "EPIC",
  plot.bias = FALSE,
  prior.prob = TRUE,
  anno = ann,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  genomic.features = "ALL",
  sig.genes = FALSE
)

df=KEGG_term[KEGG_term$P.DE<=0.05,]
write.table(KEGG_term, paste0(output_dir,"/pathway_KEGG.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# plot
# Filtering
TFTerm_KEGG=df[grepl("tuberculosis", df$TERM, ignore.case=T) |
                               grepl("interleukin", df$TERM, ignore.case=T) |
                               grepl("chemokine", df$TERM, ignore.case=T) |
                               grepl("immune", df$TERM, ignore.case=T) | 
                               grepl("kappa", df$TERM, ignore.case=T) |
                               grepl("TNF", df$TERM, ignore.case=T) |
                               grepl("inflamatory", df$TERM, ignore.case=T),]

ggplot(df,
       aes(x=reorder(Description, DE, decreasing = FALSE), y=DE)) +
  geom_bar(aes(fill = P.DE), stat="identity") +
  scale_fill_gradientn(colours = rainbow(5)) +
  labs(x = '', y = "#of genes that are differentially methylated in the term", fill = "P.DE") +
  coord_flip() +
  theme(text = element_text(size = 20),
        axis.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold', size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(output_dir,"/region_pathway_KEGG.pdf"), width = 20, height = 10)

