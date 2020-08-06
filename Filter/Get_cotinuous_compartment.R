#!/usr/bin/Rscript
# Mask intra compartment

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="Compartment definition file"),
  make_option(c("-o", "--out"), default="NA", help="output file")
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

domain_id <- rep(0, nrow(COMP))
for(n in 2:nrow(COMP)){
  domain_id[n] <- ifelse(COMP[n-1,"compartment"] == COMP[n,"compartment"] & COMP[n-1,"chr"] == COMP[n,"chr"] , domain_id[n-1], domain_id[n-1]+1)
}

COMP <- data.frame(COMP, domain_id)

D_out <- COMP %>% dplyr::group_by(domain_id, chr,compartment) %>% 
  dplyr::summarise(start=min(start), end=max(end))

options(scipen=10)
FILE_out <- as.character(opt["out"])
write.table(D_out %>% select(chr, start, end, compartment, domain_id), FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





