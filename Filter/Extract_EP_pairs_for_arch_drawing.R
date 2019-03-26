#!/usr/bin/Rscript
# SignificantなピークからE-Pの組み合わせを出力する


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="peak file"),
  make_option(c("-o", "--out"), default="NA", help="output file"),
  make_option(c("-p", "--promoter"), default="NA", help="active promoter file"),
  make_option(c("-e", "--enhancer"), default="NA", help="enhancer file")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_sig <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
FILE_promoter <- as.character(opt["promoter"])
FILE_enhancer <- as.character(opt["enhancer"])

suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(dplyr)))


### pattern
pattern_check <- function(v){
  e <- v[1]
  p <- v[2]
  if(e + p == 0){
    pp <- "N"
  }else if(e + p == 2){
    pp <- "P"
  }else{
    pp <- ''
    if(e==1){
      pp <- paste0(pp, "E")
    }
    if(p==1){
      pp <- paste0(pp, "P")
    }
  }
  pp
}

D_out <- c()
D_sig <- read.table(FILE_sig, header=TRUE, sep="\t", stringsAsFactors = FALSE)
G_L <- GRanges(D_sig %>% dplyr::rename(chr=chr1, start=start1, end=end1) %>% select(chr, start, end))
G_R <- GRanges(D_sig %>% dplyr::rename(chr=chr2, start=start2, end=end2) %>% select(chr, start, end))

getGrange_sig <- function(file, name){
  df <- read.table(file, header=FALSE, sep="\t", stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "name", "score", "strand"))
  gg <- GRanges(df %>% select(chr, start, end))
  f1 <- countOverlaps(G_L, gg)
  f2 <- countOverlaps(G_R, gg)
  f1 <- ifelse(f1 > 0, 1, 0)
  f2 <- ifelse(f2 > 0, 1, 0)
  df <- data.frame(L=f1, R=f2, stringsAsFactors = FALSE)
  colnames(df) <- paste0(name, "_", c("L", "R"))
  rm(gg, f1, f2)
  df
}

D_table2 <- getGrange_sig(FILE_promoter, "E")
D_table2 <- cbind(D_table2, getGrange_sig(FILE_enhancer, "P"))

LLL <- apply(D_table2[, c("E_L", "P_L")], 1, pattern_check)
RRR <- apply(D_table2[, c("E_R", "P_R")], 1, pattern_check)
rm(D_table2)
category <- c("E", "P","N")
D_table <- data.frame(LLL=factor(LLL, levels=category), RRR=factor(RRR, levels = category))
rm(LLL, RRR)
D_table <- D_table %>% mutate(comb=if_else(as.integer(LLL) < as.integer(RRR), paste(LLL, RRR, sep="-"), paste(RRR, LLL, sep="-")))


comb_order <- unique(D_table$comb)
comb_order <- comb_order[-which(comb_order=="E-P")]
comb_order <- c(comb_order, "E-P")

D_out <- data.frame(D_sig, comb=factor(D_table$comb, levels = comb_order)) %>% arrange(comb)
write.table(D_out, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)









