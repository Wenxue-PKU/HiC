#!/usr/bin/Rscript
# Convert HOMER matrices to original format

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="HOMER output matrices"),
  make_option(c("-o", "--out"), default="NA", help="output rds file")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])

# FILE_in <- "X:/hideki_projects/388_20180709_HiC_Tempera_Human/data2/LMP1-con/40kb/ICE/chr20.txt"


DATA_head <- read.table(FILE_in, header=TRUE, sep="\t", nrows = 5)
classes <- sapply(DATA_head, class)
classes[1] <- "NULL"
map <- read.table(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE, colClasses = classes, row.names = "Regions")
rm(DATA_head, classes)


r <- rownames(map)
LocList <- strsplit(r, "-")
LocMatrix <- matrix(unlist(LocList), ncol=2, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[2,2]) - as.numeric(LocMatrix[1,2])
options(scipen=10)
r <- paste(LocMatrix[,1], LocMatrix[,2], as.integer(LocMatrix[,2]) + Resolution -1, sep=":")

colnames(map) <- r
rownames(map) <- r
map <- as.matrix(map)
map[is.infinite(map)] <- NA

saveRDS(map, FILE_out)



