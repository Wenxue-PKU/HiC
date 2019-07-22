#!/usr/bin/Rscript
# Circos plot


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="loop list"),
  make_option(c("-o", "--out"), default="NA", help="output image"),
  make_option(c("-g", "--gene"), default="NA", help="gene information")
)
opt <- parse_args(OptionParser(option_list=option_list))



suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(migest)))
suppressWarnings(suppressMessages(library(dplyr)))


# FILE_chromosome <- as.character(opt["length"])
# FILE_chromosome <- "W:/Data/Human_seq/hg19_EBV/LENGTH.txt"
FILE_chromosome <- "/wistar/noma/Data/Human_seq/hg19_EBV/LENGTH.txt"
chromosomes <- c("EBV")

FILE_gene <- as.character(opt["gene"])


data_chr <- list()
data_chr_tmp <- read.table(FILE_chromosome, header=FALSE, sep="\t", stringsAsFactors = FALSE)
rownames(data_chr_tmp) <- data_chr_tmp[,1]
data_chr_tmp <- data_chr_tmp[chromosomes,]
data_chr$df <- data.frame(chr=factor(data_chr_tmp[,1], levels=chromosomes), start=0, end=data_chr_tmp[,2])
rm(data_chr_tmp)


D_gene <- read.table(FILE_gene, header=TRUE, sep="\t", stringsAsFactors = FALSE)
D_gene <- D_gene %>% mutate(chr="EBV") %>% select(chr, start, end, Symbol)

FILE_in <- as.character(opt["in"])
df <- read.table(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE)
df <- df %>% filter(qval < 0.05)
data_link <- list()
data_link$bed1 <- df %>% select(start1, end1) %>% rename(start=start1, end=end1) %>% mutate(chr="EBV", start=start + 2000, end=end-2000) %>% select(chr, start, end)
data_link$bed2 <- df %>% select(start2, end2) %>% rename(start=start2, end=end2) %>% mutate(chr="EBV", start=start + 2000, end=end-2000) %>% select(chr, start, end)



circos.clear()

FILE_out <- as.character(opt["out"])
pdf(file =FILE_out, height=10, width=10)
circos.par(track.height = 0.1, cell.padding = c(0,0,0,0), track.margin = c(0,0), start.degree = 90, gap.degree =4)
circos.genomicInitialize(data_chr$df, axis.labels.cex = 1)
circos.genomicLabels(D_gene, labels.column = 4, side="outside")
circos.genomicTrack(data_chr$df, ylim=c(0,1),
                    panel.fun = function(region, value, ...){
                      
                      circos.genomicRect(D_gene %>% select(start, end), ytop=1, ybottom=0, col=adjustcolor('orange', alpha.f = 0.3))
                    })

circos.genomicLink(data_link$bed1, data_link$bed2, col = adjustcolor('blue', alpha.f = 0.2))
circos.clear()
dummy <- dev.off()






