#!/usr/bin/Rscript
# Extract input file location's score from HiC matrices

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="pair information to extract (chr1,start1,end1,chr2,start2,end2)"),
  make_option(c("--dir"), default="NA", help="directory of matrices files"),
  make_option(c("--normalize"), default=FALSE, help="distance normalize (TURE) or not"),
  make_option(c("-o", "--out"),help="output files of scores"),
  make_option(c("--header"), default=FALSE, help="ouput location information or not")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_in <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2018-10-02_extract_significantly_different_pairs/data2_OISx1_40kb/top50_pairs.txt"
DIR <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/data2/IMR90_OIS_con/40kb/Raw/"

options(scipen=10)
FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
DIR <- as.character(opt["dir"])
if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
  paste(DIR, "/", sep="")
}


suppressWarnings(suppressMessages(library(dplyr)))

D_target <- read.table(FILE_in, header=FALSE, sep="\t", stringsAsFactors = FALSE)
D_target <- D_target[,1:6]
colnames(D_target) <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
D_target <- D_target %>% mutate(chr=chr1, key1=paste(chr1, start1, end1, sep=":"), key2=paste(chr2, start2, end2, sep=":")) %>% select(chr, key1, key2)

D_score <- c()
Chromosomes <- unique(D_target %>% pull(chr))
# options(check.bounds = FALSE)
for(c in Chromosomes){
  FILE_object <- paste(DIR, c, ".rds", sep="")
  if(file.exists(FILE_object)){
    map <- readRDS(FILE_object)
    r <- rownames(map)
    # LocList <- strsplit(r, ":")
    # LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
    # RESOLUTION <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
    
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
    
    D_score_c <- D_target %>% filter(chr==c & key1 %in% r & key2 %in% r) %>% select(-chr)
    
    if(nrow(D_score_c) > 0){
      score <- map[D_score_c %>% as.matrix()]
      D_score_c <- D_score_c %>% mutate(score=score)
      if(is.null(D_score)){
        D_score <- D_score_c
      }else{
        D_score <- rbind(D_score, D_score_c)
      }
    }
  }
}
D_target <- dplyr::left_join(D_target, D_score, by=c("key1", "key2"))
if(as.character(opt["header"]) == "TRUE"){
  write.table(D_target %>% select(key1, key2, score), FILE_out, row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
}else{
  write.table(D_target %>% select(score), FILE_out, row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
}

