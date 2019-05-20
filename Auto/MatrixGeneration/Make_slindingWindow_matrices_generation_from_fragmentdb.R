#!/usr/bin/Rscript
# Create intra-chromosome sliding window matrices from database

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--db"), default="NA", help="fragment db"),
  make_option(c("-r", "--resolution"), default="NA", help="resolution (bp)"),
  make_option(c("-s", "--sliding"), default="NA", help="sliding window (bp)"),
  make_option(c("-o", "--out"), default="NA", help="matrices name xxx.matrix"),
  make_option(c("-c", "--chr"), default="NA", help="target chromosome"),
  make_option(c("--start"), default="NA", help="start area"),
  make_option(c("--end"), default="NA", help="end area")
)
opt <- parse_args(OptionParser(option_list=option_list))

options(scipen=10)
suppressWarnings(suppressMessages(library(dplyr)))


#=====================================
# get path of program directory
#=====================================
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
PROGRAM_getmatrix <- paste0(script.basename, "/Make_association_from_fragmentdb_selectedArea.pl")


DB_fragment <- as.character(opt["db"])
RESOLUTION <- as.numeric(as.character(opt["resolution"]))
SLIDING <- as.character(opt["sliding"])
if(SLIDING == "NA"){
  SLIDING <- RESOLUTION
}else{
  SLIDING <- as.numeric(SLIDING)
}
SURROUNDING_BIN <- ((RESOLUTION / SLIDING) - 1)/2

TARGET_CHR <- as.character(opt["chr"])
TARGET_START <- as.character(opt["start"])
TARGET_END <- as.character(opt["end"])

FILE_matrix <- as.character(opt["out"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
FILE_tmp <- paste0(FILE_matrix, "_tmp_matrixcalc")
system(paste("perl", PROGRAM_getmatrix, "-i", DB_fragment, "-o", FILE_tmp, "-r", SLIDING, "-c", TARGET_CHR, "-s", TARGET_START, "-e", TARGET_END))

map <- as.matrix(read.table(FILE_tmp, header=TRUE, check.names = FALSE))

NUM_LINE <- nrow(map)
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Location <- data.frame(id=1:NUM_LINE, chr=LocMatrix[,1], start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]), stringsAsFactors = FALSE)
rm(r, LocList, LocMatrix)


SURROUNDING_BIN <- ((RESOLUTION / SLIDING) - 1)/2

if(SURROUNDING_BIN != 0){
  map_new <- map
  for (i in (1:nrow(map))){
    CalcAverage <- function(m){
      sum(map[max(1,i-SURROUNDING_BIN):min(nrow(map),i+SURROUNDING_BIN),
              max(1,m-SURROUNDING_BIN):min(ncol(map),m+SURROUNDING_BIN)], na.rm = TRUE)
    }
    
    score <- sapply(i:nrow(map), CalcAverage)
    map_new[cbind(i, i:nrow(map))] <- score
    map_new[cbind(i:nrow(map), i)] <- score
  }
  map <- map_new
  rm(map_new)
}

Location <- Location %>% mutate(ss=ifelse(start - SURROUNDING_BIN * SLIDING < Location[1,"start"], Location[1,"start"], start - SURROUNDING_BIN * SLIDING),
                                ee=ifelse(end + SURROUNDING_BIN * SLIDING > Location[nrow(Location),"end"], Location[nrow(Location),"end"], end + SURROUNDING_BIN * SLIDING))
Location <- Location %>% mutate(binname=paste(chr, ss, ee, sep=":")) %>% pull(binname)
colnames(map) <- Location
rownames(map) <- Location
rm(Location)


saveRDS(map, FILE_object)
write.table(map, FILE_matrix, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

file.remove(FILE_tmp)




