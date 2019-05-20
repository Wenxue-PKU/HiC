#!/usr/bin/Rscript
# Extract total scores for specific distance ranges

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="hic matrix"),
  make_option(c("--min"), default="0", help="minimum distance"),
  make_option(c("--max"), default="999999999", help="maximum distance"),
  make_option(c("--inter"), default="FALSE", help="inter-chromosome score is add ? or not")
)
opt <- parse_args(OptionParser(option_list=option_list))



Threshold_min <- as.numeric(as.character(opt["min"]))
Threshold_max <- as.numeric(as.character(opt["max"]))


FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}


Total_score <- 0

r <- rownames(map)
LocList <- strsplit(r, ":")
if(length(LocList[[1]]) == 3){
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
  NUM_LINE <- nrow(map)
  
  for(d in 0:(NUM_LINE-1)){
    distance <- d*Resolution
    if(distance >= Threshold_min & distance <= Threshold_max){
      index1 <- 1:(NUM_LINE - d)
      index2 <- index1 + d
      index3 <- cbind(index1, index2)
      index_intra <- as.character(LocMatrix[index1,1]) == as.character(LocMatrix[index2,1])
      index3 <- index3[index_intra,]
      total <- sum(as.numeric(map[index3]), na.rm = TRUE)
      Total_score <- Total_score + total
    }
  }
}else{
  LocMatrix <- data.frame(chromosome=r)
}



if(eval(as.character(opt["inter"]))){
  CHRs <- unique(as.character(LocMatrix[,1]))
  comb <- t(combn(CHRs, 2, simplify = TRUE))
  for(i in 1:nrow(comb)){
    index1 <- which((as.character(LocMatrix[,1]) == as.character(comb[i,1])))
    index2 <- which((as.character(LocMatrix[,1]) == as.character(comb[i,2])))
    index3 <- cbind(index1, index2)
    total <- sum(as.numeric(map[index3]), na.rm = TRUE)
    Total_score <- Total_score + total
  }
}


cat(Total_score, "\n")


