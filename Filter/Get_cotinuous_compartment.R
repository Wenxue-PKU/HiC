#!/usr/bin/Rscript
# Mask intra compartment

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="Compartment definition file"),
  make_option(c("-o", "--out"), default="NA", help="mask object")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(dplyr))

FILE_comp <- as.character(opt["in"])
COMP_head <- read.table(FILE_comp, header=FALSE, nrows = 5, sep="\t")
classes <- sapply(COMP_head, class)
if(length(classes) == 5){
  classes <- c("character", "integer", "integer", "NULL", "character")
}else{
  classes <- c("character", "integer", "integer", "NULL", "NULL", "character")
}

COMP <- read.table(FILE_comp, header=FALSE, sep="\t", stringsAsFactors = FALSE, colClasses = classes)
colnames(COMP) <- c("chr", "start", "end", "compartment")
COMP[is.na(COMP[,"compartment"]), "compartment"] <- "N"

index_comp <- rep(0, nrow(COMP))
for(n in 2:nrow(COMP)){
  index_comp[n] <- ifelse(COMP[n-1,"compartment"] == COMP[n,"compartment"] & COMP[n-1,"chr"] == COMP[n,"chr"] , index_comp[n-1], index_comp[n-1]+1)
}

COMP <- data.frame(COMP, index_comp)

data <- COMP %>% dplyr::group_by(index_comp, chr,compartment) %>% 
  dplyr::summarise(start=min(start), end=max(end))

options(scipen=10)
Target <- list()

Target$chr <- unique(COMP[,"chr"])
Target$A <- as.data.frame(data %>% dplyr::filter(compartment=="A") %>% ungroup() %>% select(chr, start, end))
Target$B <- as.data.frame(data %>% dplyr::filter(compartment=="B") %>% ungroup() %>% select(chr, start, end))

FILE_out <- as.character(opt["out"])
saveRDS(Target, FILE_out)


