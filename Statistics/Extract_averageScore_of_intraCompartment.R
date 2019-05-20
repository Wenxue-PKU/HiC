#!/usr/bin/Rscript
# Extract contact probabilities of specified intra-compartment

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--dir"), default="NA", help="directory of matrices"),
  make_option(c("-i", "--in"), default="NA", help="target location .rds"),
  make_option(c("-o", "--out"), default="NA", help="output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

# DIR_data <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/data/IMR90_OIS_con/200kb/ICE"
# FILE_comp <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2018-07-31_compare_contact_probabilities_at_compartment_between_samples/Region_intraCompartment.rds"


FILE_out <- as.character(opt["out"])

DIR_data <- as.character(opt["dir"])
if(substring(DIR_data, nchar(DIR_data), nchar(DIR_data)) != "/"){
  DIR_data <- paste(DIR_data, "/", sep="")
}
FILE_comp <- as.character(opt["in"])
TARGET <- readRDS(FILE_comp)
options(scipen=10)
TABLE <- c()
for(chr in TARGET$chr){
  FILE_mat <- paste(DIR_data, chr, ".rds", sep="")
  map <- readRDS(FILE_mat)
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
    map_expect[index3] <- Average
    map_expect[index4] <- Average
  }
  map <- ifelse(map_expect == 0, NA, log2(map / map_expect))
  
  
  
  for(n in which(TARGET$A[,"chr"] == chr)){
    bins <- seq(TARGET$A[n,"start"], TARGET$A[n,"end"], by = RESOLUTION)
    binNames <- paste(chr, bins, bins+RESOLUTION-1, sep=":")
    binNames <- binNames[binNames %in% r]
    scores <- map[as.matrix(expand.grid(binNames, binNames))]
    new_line <- c(chr, TARGET$A[n,"start"], TARGET$A[n,"end"], "A", mean(scores, na.rm = TRUE))
    TABLE <- rbind(TABLE, new_line)
  }
  
  for(n in which(TARGET$B[,"chr"] == chr)){
    bins <- seq(TARGET$B[n,"start"], TARGET$B[n,"end"], by = RESOLUTION)
    binNames <- paste(chr, bins, bins+RESOLUTION-1, sep=":")
    binNames <- binNames[binNames %in% r]
    scores <- map[as.matrix(expand.grid(binNames, binNames))]
    new_line <- c(chr, TARGET$B[n,"start"], TARGET$B[n,"end"], "B", mean(scores, na.rm = TRUE))
    TABLE <- rbind(TABLE, new_line)
  }
}
colnames(TABLE) <- c("chr", "start", "end", "compartment", "ave_score")
TABLE <- TABLE[order(TABLE[,"chr"], TABLE[,"start"], decreasing = FALSE), ]
write.table(TABLE, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



