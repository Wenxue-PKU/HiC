#!/usr/bin/Rscript
# check overlapped of interactions


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--file1"), default="NA", help="interaction information of file1"),
  make_option(c("-b", "--file2"), default="NA", help="interaction information of file2"),
  make_option(c("--name1"), default="", help="name of file1"),
  make_option(c("--name2"), default="", help="name of file2"),
  make_option(c("-o", "--out"), default="NA", help="output image file")
  
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(InteractionSet)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(VennDiagram)))
suppressWarnings(suppressMessages(library(Cairo)))
suppressWarnings(suppressMessages(library(grDevices)))

FILE_hic1 <- as.character(opt["file1"])
FILE_hic2 <- as.character(opt["file2"])
FILE_out <- as.character(opt["out"])
NAME_1 <- as.character(opt["name1"])
NAME_2 <- as.character(opt["name2"])


getInteraction <- function(file){
  D_hic <- read.table(file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  G_L <- GRanges(D_hic %>% dplyr::rename(chr=chr1, start=start1, end=end1) %>% select(chr, start, end))
  G_R <- GRanges(D_hic %>% dplyr::rename(chr=chr2, start=start2, end=end2) %>% select(chr, start, end))
  GI <- GInteractions(G_L, G_R)
  GI
}

Gi_1 <- getInteraction(FILE_hic1)
Gi_2 <- getInteraction(FILE_hic2)


Common <- countOverlaps(Gi_1, Gi_2, use.region="both")
Common <- sum(Common>0)
UniqA <- length(Gi_1) - Common
UniqB <- length(Gi_2) - Common


p <- venn.diagram(list(A = 1:(UniqA+Common), B = (UniqA+1):(UniqA+Common+UniqB)),
             alpha = 0.5, cex = 3, cat.fontface = "bold",lty =2, cat.cex = 2, fontface = "bold", fill = c("#91cf60", "#fc8d59"),
             filename = NULL, category.names = c(NAME_1, NAME_2), margin=0.02, imagetype = 'png')
CairoPNG(FILE_out, width=15, height=15, units = "cm", res=200)
grid.draw(p)
d <- dev.off()




