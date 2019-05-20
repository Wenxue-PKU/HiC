#!/usr/bin/Rscript
# BEDファイルからPairリストとOverlapするものを抽出する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--BED"), default="NA", help="BED file"),
  make_option(c("-p", "--pair"), default="NA", help="pair file. at least chr1, start1, end1, chr2, start2, end2 are required"),
  make_option(c("-o", "--out"), default="NA", help="Filtered list or add additional column with column name"),
  make_option(c("-n", "--name"), default="Overlap", help="if add column indicating Overlap(1) or not(0), name of column. NA if not output"),
  make_option(c("--header"), default="FALSE", help="output header or not"),
  make_option(c("--onlyOverlap"), default="FALSE", help="output only overlapped BED")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_bed <- as.character(opt["BED"])
FILE_pair <- as.character(opt["pair"])
FILE_out <- as.character(opt["out"])
COLUMN_NAME <- as.character(opt["name"])
FLAG_onlyOverlap <- as.character(opt["onlyOverlap"])
FLAG_header <- eval(parse(text=as.character(opt["header"])))

suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(dplyr)))


D_pair <- read.table(FILE_pair, header=TRUE, sep="\t", stringsAsFactors = FALSE)
G_L <- GRanges(D_pair %>% dplyr::rename(chr=chr1, start=start1, end=end1) %>% select(chr, start, end))
G_R <- GRanges(D_pair %>% dplyr::rename(chr=chr2, start=start2, end=end2) %>% select(chr, start, end))

D_bed <- read.table(FILE_bed, header=FALSE, sep="\t", stringsAsFactors = FALSE)
colnames(D_bed)[1:3] <-  c("chr", "start", "end")
gg <- GRanges(D_bed %>% select(chr, start, end))
f1 <- countOverlaps(gg, G_L)
f2 <- countOverlaps(gg, G_R)
score <- ifelse(f1+f2 > 0, 1, 0)


D_bed <- D_bed %>% mutate(score=score)

if(FLAG_onlyOverlap == "TRUE"){
  D_bed <- D_bed %>% filter(score == 1)
}

if(COLUMN_NAME == "NA"){
  D_out <- D_bed %>% select(-score)
}else{
  D_out <- D_bed %>% dplyr::rename(!!COLUMN_NAME := score)
}

write.table(D_out, FILE_out, sep = "\t", row.names = FALSE, col.names = FLAG_header, quote = FALSE)






