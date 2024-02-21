# Check installed packages
if (!require(c("FactoMineR", "factoextra","ggplot2"))) {
  install.packages(c("FactoMineR", "factoextra","ggplot2"))
}

# Load libraries
library("FactoMineR")
library("factoextra")

# Load fixed dataset
cesc <- read.table('database/fixed/fixed-gdc-tcga-cesc.csv', header = TRUE, sep = ",",na.strings=c(NA,''))
hnsc <- read.table('database/fixed/fixed-gdc-tcga-hnsc.csv', header = TRUE, sep = ",",na.strings=c(NA,''))
lihc <- read.table('database/fixed/fixed-gdc-tcga-lihc.csv', header = TRUE, sep = ",",na.strings=c(NA,''))