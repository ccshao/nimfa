## File Name : nmfHeatmap.R
## Created By : Chunxuan Shao
## Creation Date : Mon Jan 26 15:01:35 2015
## Last Modified : Fri 26 May 2017 05:11:29 PM CEST
## Description: 1.) draw heatmap from cell basis matrix;

library(RColorBrewer)
## NMF does not work now in R 3.4.0
suppressMessages(library(vegan))
library(pheatmap)

################################################################################
rsplit <- function (x, s, n) {
    p <- paste0("[^", s, "]*")
        rx <- paste0(s, "(?=", paste(rep(paste0(p, s), n - 1), collapse = ""), 
                        p, "$)")
    unlist(strsplit(x, rx, perl = TRUE))
}

ff_basisMatrix <- function(fileName) {
  ## read the normalized basis matrix;
  ## plot the heatmap;
  # browser()
  dat.nor <- read.delim(fileName, row.names = 1)
  filePre <- rsplit(fileName, ".", 1)[1]
  # cat("+++ Making heatmap for", filePre, "+++\n")
  
  pdf(paste(filePre, ".metaCells.heatmap.pdf", sep = ""), w = 5, h = 12, onefile=FALSE)
  # xx.1 <- aheatmap_wrapper(dat.nor, col = colorRampPalette(rev(brewer.pal(11, "PuOr")))(100), fontsize = 12)
  # xx.1 <- aheatmap(dat.nor, col = colorRampPalette(rev(brewer.pal(11, "PuOr")))(100), fontsize = 12)
  ## to reorder the dendgram as heatmap in base
  xx.1 <- pheatmap(dat.nor, 
                   col = colorRampPalette(rev(brewer.pal(11, "PuOr")))(100), 
                  show_rownames = FALSE )
  dev.off()
  write.table(rownames(dat.nor)[xx.1$tree_row$order], file = paste0(filePre, ".cells.order.inheatmap.csv"), sep = "\t", quote = FALSE)
  # xx.2 <- heatmap(dat.nor, scale = "none")
  # write.table( rev(rownames(dat.nor)[xx.2$rowInd]), file = "1.csv", sep = "\t", quote = FALSE)

}

################################################################################

args <- commandArgs(trailingOnly = TRUE)
# print(args)
fileName <- args[1]

ff_basisMatrix(fileName)

