#!/usr/bin/Rscript
# Extract input file location's score from HiC matrices

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="pair information to extract (chr1,start1,end1,chr2,start2,end2)"),
  make_option(c("--dir"), default="NA", help="directory of matrices files"),
  make_option(c("-o", "--out"),help="output files of scores"),
  make_option(c("--na.rm"), default=FALSE, help="remove NA from output"),
  make_option(c("--force_output"), default=FALSE, help="even output file exists proceed")
)
opt <- parse_args(OptionParser(option_list=option_list))


# FILE_in <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2018-08-31_association_between_protein_binding_sites/combinations/CAPH2_275_S2_200kb_0_500kb.txt"
# DIR <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/data/IMR90_OIS_con/200kb/Raw/"

options(scipen=10)
FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
DIR <- as.character(opt["dir"])
if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
  paste(DIR, "/", sep="")
}

### もし出力ファイルが既に存在していたら、中止
if(as.character(opt["force_output"]) != "TRUE"){
  if(file.exists(FILE_out)){
    cat("Output file is already exists\n")
    q()
  }
}


DATA <- read.table(FILE_in, header=FALSE, sep="\t", stringsAsFactors = FALSE)

Chromosomes <- unique(DATA[,1])
target <- list()
for(c in Chromosomes){
  index <- DATA[,1] == c
  B1 <- paste(DATA[index,1], DATA[index,2], DATA[index,3], sep=":")
  B2 <- paste(DATA[index,4], DATA[index,5], DATA[index,6], sep=":")
  target[[c]] <- cbind(B1, B2)
}
rm(DATA, B1, B2, c)

options(check.bounds = FALSE)
for(c in Chromosomes){
  FILE_object <- paste(DIR, c, ".rds", sep="")
  if(file.exists(FILE_object)){
    map <- readRDS(FILE_object)
    r <- rownames(map)
    LocList <- strsplit(r, ":")
    LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
    RESOLUTION <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
    
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
    
    ok_pair <- (target[[c]][,1] %in% rownames(map)) & (target[[c]][,2] %in% colnames(map))
    score <- rep(NA, nrow(target[[c]]))
    score[ok_pair] <- map[target[[c]][ok_pair,]]
    
    ### NA を除く（ただし、この場合場所指定ファイルとの整合性が取れないことに注意)
    if(as.character(opt["na.rm"]) == "TRUE"){
      score <- score[!is.na(score)]
    }
    write.table(score, FILE_out, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  }
}


