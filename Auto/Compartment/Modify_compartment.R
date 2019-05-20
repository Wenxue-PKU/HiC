#!/usr/bin/Rscript
# compartmentの修正プログラム

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="original PCA file"),
  make_option(c("-o", "--out"), help="modified file")
)
opt <- parse_args(OptionParser(option_list=option_list))


DATA <- read.table(as.character(opt["in"]), sep="\t", stringsAsFactors = FALSE, header=FALSE, colClasses = c("character", "numeric", "numeric", "numeric","character"), na.strings = "")

Edit_Compartment <- function(Com){
  Com.new <- Com
  Previous_status <- "NA"
  Current_length <- 1
  Current_status <- Com[1]
  Current_start <- 1
  for(i in 2:length(Com)){
    if(Com[i] == Current_status){
      Current_length <- Current_length + 1
    }else{
      # binの長さが5bin未満であれば、前のstatusを継続する
      if(Current_length < 5){
        Com.new[Current_start:(i-1)] <- Previous_status
      }else{
        Previous_status <- Current_status
      }
      Current_length <- 1
      Current_start <- i
      Current_status <- Com[i]
    }
  }
  Com.new
}

Compartment <- as.character(DATA[,5])
Compartment.new <- Compartment

Chromosomes <- unique(as.character(DATA[,1]))
for(c in Chromosomes){
  index_chr <- as.character(DATA[,1]) == c
  Compartment.new[index_chr] <- Edit_Compartment(Compartment[index_chr])
}

options(scipen=10)
write.table(cbind(DATA, Compartment.new), file=as.character(opt["out"]), sep="\t", eol = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)





