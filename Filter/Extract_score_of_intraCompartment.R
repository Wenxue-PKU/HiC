#!/usr/bin/Rscript
# Extract contact probabilities of specified intra-compartment

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--dir"), default="NA", help="directory of matrices"),
  make_option(c("-i", "--in"), default="NA", help="target location .rds"),
  make_option(c("-o", "--out"), default="NA", help="prefix of output text file(output is **_A.txt, **_B.txt)")
)
opt <- parse_args(OptionParser(option_list=option_list))

# DIR_data <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/data/IMR90_OIS_con/200kb/ICE"
FILE_comp <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2018-07-31_compare_contact_probabilities_at_compartment_between_samples/Region_intraCompartment.rds"


FILE_out_prefix <- as.character(opt["out"])
FILE_out_A <- paste(FILE_out_prefix, "_A.txt", sep="")
FILE_out_B <- paste(FILE_out_prefix, "_B.txt", sep="")
if(file.exists(FILE_out_A)){
  file.remove(FILE_out_A)
}
if(file.exists(FILE_out_B)){
  file.remove(FILE_out_B)
}

DIR_data <- as.character(opt["dir"])
if(substring(DIR_data, nchar(DIR_data), nchar(DIR_data)) != "/"){
  DIR_data <- paste(DIR_data, "/", sep="")
}
FILE_comp <- as.character(opt["in"])
TARGET <- readRDS(FILE_comp)


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
    scores <- scores[!is.na(scores)]
    write.table(scores, FILE_out_A, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  }
  
  for(n in which(TARGET$B[,"chr"] == chr)){
    bins <- seq(TARGET$B[n,"start"], TARGET$B[n,"end"], by = RESOLUTION)
    binNames <- paste(chr, bins, bins+RESOLUTION-1, sep=":")
    binNames <- binNames[binNames %in% r]
    scores <- map[as.matrix(expand.grid(binNames, binNames))]
    scores <- scores[!is.na(scores)]
    write.table(scores, FILE_out_B, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
  }
}


 

