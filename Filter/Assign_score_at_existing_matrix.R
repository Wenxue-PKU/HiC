#!/usr/bin/Rscript
# 既存のmatrixをベースに与えたpairのスコアを代入したマトリックスを出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="score list for assign (key1, key2, score)"),
  make_option(c("--temp"), default="NA", help="template matrix"),
  make_option(c("-o", "--out"),help="output matrix")
)
opt <- parse_args(OptionParser(option_list=option_list))

options(scipen=10)


FILE_matrix <- as.character(opt["temp"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(file.exists(FILE_object)){
  map <- readRDS(FILE_object)
}else if(file.exists(FILE_matrix)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  cat(FILE_matrix, " was not found\n")
  q()
}
map[,] <- NA


D_target <- read.table(as.character(opt["in"]), header=FALSE, sep="\t", stringsAsFactors = FALSE, check.names = F)
nc <- colnames(map)
nr <- rownames(map)
OK_pair <- D_target[,1] %in% nr & D_target[,2] %in% nc
D_target <- D_target[OK_pair,]


target <- cbind(D_target[,1], D_target[,2])
map[target] <- D_target[,3]


write.table(map, file=as.character(opt["out"]), quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA)




