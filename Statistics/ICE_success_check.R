#!/usr/bin/Rscript
# ICE normalization success check. Count % of NA lines

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("--method"), default="p", help="p(percent) or n(number)")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_matrix <- as.character(opt["in"])
FILE_object <- sub("matrix", "rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

if(as.character(opt["method"]) == "p"){
	cat(sum(apply(map, 1, sum, na.rm=TRUE) == 0) / nrow(map), "\n")
}else{
	cat(sum(apply(map, 1, sum, na.rm=TRUE) == 0), "\n")
}

