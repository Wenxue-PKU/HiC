#!/usr/bin/Rscript
# 2つのBEDファイルで指定された領域のすべての組み合わせ領域のHi-Cスコアを取得する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--data"),help="matrix directory"),
  make_option(c("--resolution"),help="matrix resolution"),
  make_option(c("-o", "--out"),help="output file"),
  make_option(c("-a", "--bed1"),help="bed file for region1"),
  make_option(c("-b", "--bed2"),help="bed file for region2"),
  make_option(c("--max_distance"), default=500000, help="maximum distance between two region")
)
opt <- parse_args(OptionParser(option_list=option_list))

DIR_data <- as.character(opt["data"])
if(substring(DIR_data, nchar(DIR_data), nchar(DIR_data)) != "/"){
  DIR_data <- paste(DIR_data, "/", sep="")
}
FILE_OUT <- as.character(opt["out"])
Resolution <- as.numeric(as.character(opt["resolution"]))


options(scipen=10)
getBEDregion <- function(file){
  DATA_head <- read.table(file, sep="\t", header=FALSE, nrows = 5)
  COLUMN_NUM <- ncol(DATA_head)
  DATA_all <- read.table(file, header=FALSE, colClasses= c("character", "numeric", "numeric", "character", rep("NULL", COLUMN_NUM - 4)), stringsAsFactors = FALSE, sep = "\t", check.names = FALSE)
  Middle <- (as.numeric(DATA_all[,2]) + as.numeric(DATA_all[,3]))/2
  Bin_start <- as.integer(Middle / Resolution) * Resolution
  Bin_end <- Bin_start + Resolution - 1
  Bin_id <- paste(DATA_all[,1], Bin_start, Bin_end, sep=":")
  BED <- cbind(DATA_all, Bin_id)
  BED
}

BED1 <- getBEDregion(as.character(opt["bed1"]))
BED2 <- getBEDregion(as.character(opt["bed2"]))
Chromosomes <- unique(as.character(BED1[,1]))

for(CHR in Chromosomes){
  map <- readRDS(paste(DIR_data, CHR, ".rds", sep=""))
  r <- rownames(map)
  
  BED1.chr <- BED1[as.character(BED1[,1]) == CHR,]
  BED2.chr <- BED2[as.character(BED2[,1]) == CHR,]
  BED1.chr[!(BED1.chr[,5] %in% r),5] <- NA
  BED2.chr[!(BED2.chr[,5] %in% r),5] <- NA

  COMB <- expand.grid(rownames(BED1.chr), rownames(BED2.chr), stringsAsFactors = FALSE)
  
  # distance restriction
  COMB <- COMB[abs(as.numeric(BED1.chr[COMB[,1],2]) - as.numeric(BED2.chr[COMB[,2],2])) < as.numeric(as.character(opt["max_distance"])),]
  
  pair <- cbind(as.character(BED1.chr[COMB[,1],5]), as.character(BED2.chr[COMB[,2],5]))
  HiC_score <- map[pair]
  
  APPEND=TRUE
  if(CHR == Chromosomes[1]){
    APPEND=FALSE
  }
  suppressWarnings(write.table(cbind(BED1.chr[COMB[,1],1:4], BED2.chr[COMB[,2],1:4], HiC_score), file = FILE_OUT,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", eol = "\n", append=APPEND))
}



