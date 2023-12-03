library(VennDiagram)
library(RColorBrewer)

data_dir="/path/to/data_dir"

myCol <- brewer.pal(4, "Set1")

set1 <- read.table(paste0(data_dir, "/data/train/DMP/TB-HC/BMA.txt"), header=T)$probe_ID[-1]
set2 <- read.table(paste0(data_dir, "/data/train/DMP/TB-HC/RF_min.txt"), header=T)$name
set3 <- read.table(paste0(data_dir, "/data/train/DMP/TB-HC/PLS_min.txt"), header=T)$name
set4 <- read.table(paste0(data_dir, "/data/train/DMP/TB-HC/LASSO.txt"), header=T)$probeID

venn.diagram(
x = list(set1, set2, set3, set4),
category.names = c(paste0("BMA\n(", length(set1), ")" ) , paste0("MUVR - RF\n(", length(set2), ")"), paste0("MUVR - PLS\n(", length(set3), ")"), paste0("LASSO\n(", length(set4), ")")),
filename = paste0(data_dir, "/data/train/DMP/TB-HC/14_venn_diagramm.pdf"),
output=TRUE,

# Circles
lwd = 3,
lty = 'blank',
fill = myCol,
col = "black",

# Numbers
cex = 1,
fontface = "bold",
fontfamily = "sans",
)

