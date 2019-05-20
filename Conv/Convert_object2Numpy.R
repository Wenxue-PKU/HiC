#!/usr/bin/Rscript
# Convert matrix object to numpy format

library(RcppCNPy)

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="rds file"),
  make_option(c("-o", "--out"), default="NA", help="output numpy file")
)
opt <- parse_args(OptionParser(option_list=option_list))


mat <- readRDS(as.character(opt["in"]), "integer")
colnames(mat) <- NULL
rownames(mat) <- NULL

index <- cbind(1:nrow(mat), 1:nrow(mat))
mat[index] <- 0
mat <- round(mat)

npySave(as.character(opt["out"]), mat, mode='w', checkPath = TRUE)
