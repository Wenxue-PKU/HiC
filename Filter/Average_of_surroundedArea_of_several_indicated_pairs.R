#!/usr/bin/Rscript
# Make average matrix of surrounded area of indicated pairs

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="pair information to extract (chr1,start1,end1,name1,chr2,start2,end2,name2)"),
  make_option(c("--dir"), default="NA", help="directory of matrices files"),
  make_option(c("--normalize"), default=FALSE, help="distance normalize (TURE) or not"),
  make_option(c("--total_adjust"), default="NA", help="adjust total read number"),
  make_option(c("-o", "--out"),help="output file of matrix"),
  make_option(c("--bin"), default="10", help="bin size to capture from center")
)
opt <- parse_args(OptionParser(option_list=option_list))

options(scipen=10)
FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
DIR <- as.character(opt["dir"])
if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
  DIR <- paste(DIR, "/", sep="")
}
BinNum <- as.numeric(as.character(opt["bin"]))

suppressWarnings(suppressMessages(library(dplyr)))

D_target <- read.table(FILE_in, header=FALSE, sep="\t", stringsAsFactors = FALSE)
D_target <- D_target[,1:8]
colnames(D_target) <- c("chr1", "start1", "end1", "name1", "chr2", "start2", "end2", "name2")
D_target <- D_target %>% mutate(chr=chr1, key1=paste(chr1, start1, end1, sep=":"), key2=paste(chr2, start2, end2, sep=":")) %>% select(chr, key1, key2, name1, name2)


MatrixSize <- 2*BinNum + 1
score <- matrix(0, nrow=MatrixSize, ncol=MatrixSize)
Total_added <- 0
Chromosomes <- unique(D_target %>% pull(chr))
for(c in Chromosomes){
  FILE_matrix <- paste0(DIR, c, ".matrix")
  FILE_object <- sub(".matrix", ".rds", FILE_matrix)
  if(file.exists(FILE_object)){
    map <- readRDS(FILE_object)
  }else if(file.exists(FILE_matrix)){
    map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  }else{
    cat(FILE_object, " was not found. Skipped...\n")
    next
  }
  r <- 1:nrow(map)
  names(r) <- rownames(map)
  
  if(as.character(opt["total_adjust"]) == "TRUE"){
    map <- map / sum(map, na.rm = TRUE) * sum(!is.na(map), na.rm=TRUE)
  }
  
  
  if(as.character(opt["normalize"]) =="TRUE"){
    total_per_line <- apply(map, 1, sum)
    map[total_per_line==0,] <- NA
    map[,total_per_line==0] <- NA
    
    # Observed / Expectのmatrixに変換する
    map_expect <- map
    NUM_LINE <- nrow(map)
    for(d in 0:(NUM_LINE-1)){
      index1 <- 1:(NUM_LINE - d)
      index2 <- index1 + d
      index3 <- cbind(index1, index2)
      index4 <- cbind(index2, index1)
      Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
      if(is.na(Average)){
        Average = NA
      }else if(Average == 0){
        Average = NA
      }
      map_expect[index3] <- Average
      map_expect[index4] <- Average
    }
    map <- ifelse(map == 0, 0, log2(map / map_expect))
  }
  
  D_score_c <- D_target %>% filter(chr==c & key1 %in% names(r) & key2 %in% names(r)) %>% select(key1, key2) %>% distinct(key1, key2, .keep_all = TRUE)
  
  D_score_c2 <- data.frame(bin1=r[D_score_c$key1], bin2=r[D_score_c$key2], stringsAsFactors = FALSE)
  D_score_c2 <- D_score_c2 %>% filter(bin1 - BinNum > 0 & bin2 - BinNum > 0 & bin1 + BinNum < length(r) & bin2 + BinNum < length(r))
  
  
  
  if(nrow(D_score_c2) > 0){
    for(i in seq(-BinNum, BinNum)){
      for(j in seq(-BinNum, BinNum)){
        LL <- cbind(D_score_c2$bin1 + i, D_score_c2$bin2 + j)
        score[i+BinNum+1,j+BinNum+1] <- score[i+BinNum+1,j+BinNum+1] + sum(map[LL], na.rm = TRUE)
      }
    }
    Total_added <- Total_added + nrow(D_score_c2)
  }
}

### Average
score <- score / Total_added


write.table(score, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


