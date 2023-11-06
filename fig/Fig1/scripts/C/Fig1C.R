library(ggplot2)
library(data.table)
library(RColorBrewer)

### Setup directory
directory = "/path/to/data"

input_dir=paste0(directory, "/fig/Fig1/data/C")
output_dir=paste0(directory, "/fig/Fig1/figure/C")

GO_term=data.frame(fread(paste0(input_dir,"/pathway_GO.tsv")))
# plot
# Filtering
df=GO_term[GO_term$P.DE<=0.05,]

TFTerm_GO=df[grepl("tuberculosis", df$TERM, ignore.case=T) |
                               grepl("interleukin", df$TERM, ignore.case=T) |
                               grepl("chemokine", df$TERM, ignore.case=T) |
                               grepl("immune", df$TERM, ignore.case=T) | 
                               grepl("kappa", df$TERM, ignore.case=T) |
                               grepl("TNF", df$TERM, ignore.case=T) |
                               grepl("inflamatory", df$TERM, ignore.case=T),]

# Molecular funtions. 
ggplot(TFTerm_GO[TFTerm_GO$ONTOLOGY=="MF",],
       aes(x=reorder(TERM, DE, decreasing = FALSE), y=DE)) +
  geom_bar(aes(fill = P.DE), stat="identity") +
  scale_fill_gradientn(colours = rainbow(5)) +
  labs(x = '', y = "#of genes that are differentially methylated in the term", fill = "P.DE") +
  coord_flip() +
  theme(text = element_text(size = 20),
        axis.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold', size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(output_dir,"/MolecularFunc.pdf"), width = 20, height = 10)

# Biological process
ggplot(TFTerm_GO[TFTerm_GO$ONTOLOGY=="BP",],
       aes(x=reorder(TERM, DE, decreasing = FALSE), y=DE)) +
  geom_bar(aes(fill = P.DE), stat="identity") +
  scale_fill_gradientn(colours = rainbow(5)) +
  labs(x = '', y = "#of genes that are differentially methylated in the term", fill = "P.DE") +
  coord_flip() +
  theme(text = element_text(size = 20),
        axis.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold', size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(paste0(output_dir,"/BiologicalProc.pdf"), width = 20, height = 10)
