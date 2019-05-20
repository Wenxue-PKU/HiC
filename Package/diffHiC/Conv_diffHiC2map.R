#!/usr/bin/Rscript
# Convert diffHiC result to R object map

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="diffHiC result"),
  make_option(c("-o", "--out"), default="NA", help="map object")
)
opt <- parse_args(OptionParser(option_list=option_list))

#=============================================
# test data
#=============================================
# FILE_in <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2018-10-02_extract_significantly_different_pairs/data2_OISx3_40kb/Significant_pairs.txt"
# D_diffHiC <- D_diffHiC %>% filter(seqnames1=="chr22")


suppressWarnings(suppressMessages(library(dplyr)))

FILE_in <- as.character(opt["in"])
D_diffHiC <- read.table(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE)
D_diffHiC <- D_diffHiC %>% mutate(key1=paste(seqnames1, start1, end1, sep=":"), key2=paste(seqnames2, start2, end2, sep=":"), stringsAsFactors = FALSE)

options(scipen=10)
bin_max <- D_diffHiC %>% select(start1, start2) %>% unlist() %>% as.numeric() %>% max()
bin_min <- D_diffHiC %>% select(start1, start2) %>% unlist() %>% as.numeric() %>% min()
resolution <- as.numeric(D_diffHiC[1, "end1"]) - as.numeric(D_diffHiC[1, "start1"]) + 1
bins_start <- seq(bin_min, bin_max, by=resolution)
bins <- paste(D_diffHiC[1, "seqnames1"], bins_start, bins_start+resolution-1, sep=":")

D_data <- D_diffHiC %>% arrange(start1, start2) %>% select(key1, key2, logFC)
D_data2 <- D_data %>% filter(key1!=key2) %>% mutate(tmp=key1) %>% mutate(key1=key2, key2=tmp) %>% select(-tmp)
D_data <- rbind(D_data, D_data2) 
rm(D_data2, D_diffHiC)

map <- matrix(NA, nrow=length(bins), ncol=length(bins))
rownames(map) <- bins
colnames(map) <- bins

map[D_data %>% select(key1, key2) %>% as.matrix()] <- D_data %>% pull(logFC)
saveRDS(map, as.character(opt["out"]))




