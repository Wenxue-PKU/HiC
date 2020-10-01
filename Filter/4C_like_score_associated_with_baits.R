#!/usr/bin/Rscript
# Extract 4C score for target loci

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="directory of matrices located or rds file"),
  make_option(c("-o", "--out"),help="output directory"),
  make_option(c("--location"), default="NA", help="file describe target loci. text or excel file. chr, start, end, name, bait_chr, bait_start, bait_end"),
  make_option(c("--sufix"), default="NA", help="sufix for output file name."),
  make_option(c("--normalize"), default="NA", help="NA, average: average will be 1, probability: score were divided by total read")
)
opt <- parse_args(OptionParser(option_list=option_list))

#=============================================
# test file
#=============================================
# DIR_in <- "T:/Project/027_20200616_KD_SNPH/data/RWPE1_mix/40kb/ICE2/"
# FILE_location <- "T:/Project/027_20200616_KD_SNPH/out/2020-09-09_combined_graphs/data/4C_location.xlsx"
# DIR_out <- "T:/Project/027_20200616_KD_SNPH/out/2020-09-09_combined_graphs/data/"


options(scipen=10)
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(GenomicRanges)))


### location file
# chr, start, end, name, bait_chr, bait_start, bait_end
FILE_location <- as.character(opt["location"])
if(substring(FILE_location, nchar(FILE_location)-3, nchar(FILE_location)) == "xlsx"){
  suppressWarnings(suppressMessages(library(xlsx)))
  D_location <- read.xlsx(FILE_location, sheetIndex=1, header=TRUE, stringsAsFactors = FALSE)
}else{
  D_location <- fread(FILE_location, header=TRUE, sep="\t", stringsAsFactors = FALSE) %>% as.data.frame()
}

if(!"chr" %in% colnames(D_location)){
  cat("chr should be specified")
  q()
}
if(!"start" %in% colnames(D_location)){
  cat("start should be specified")
  q()
}
if(!"end" %in% colnames(D_location)){
  cat("end should be specified")
  q()
}
if(!"bait_chr" %in% colnames(D_location)){
  cat("bait_chr should be specified")
  q()
}
if(!"bait_start" %in% colnames(D_location)){
  cat("bait_start should be specified")
  q()
}
if(!"bait_end" %in% colnames(D_location)){
  cat("bait_end should be specified")
  q()
}

### もしnameがなかったら上から順にidを付ける
if(!"name" %in% colnames(D_location)){
  D_location <- D_location %>% mutate(name=row_number())
}

### name sufix
NAME_SUFIX <- as.character(opt["sufix"])
if(NAME_SUFIX != "NA"){
  D_location$name = paste0(D_location$name, NAME_SUFIX)
}


checkDIRpath <- function(DIR){
  if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
    DIR <- paste0(DIR, "/")
  }
  DIR
}
DIR_in <- as.character(opt["in"])
DIR_out <- checkDIRpath(as.character(opt["out"]))


for(cc in D_location %>% distinct(chr) %>% pull(chr) %>% as.character()){
  if(substring(DIR_in, nchar(DIR_in)-3, nchar(DIR_in)) != ".rds"){
    DIR_in <- checkDIRpath(DIR_in)
    map <- readRDS(paste0(DIR_in, cc, ".rds"))
  }else{
    map <- readRDS(DIR_in)
  }
  
  map <- ifelse(is.infinite(map), NA, map)
  r <- rownames(map)
  LocMatrix <- as.data.frame(r) %>% mutate(name=r) %>% tidyr::separate(r, c("chr", "start", "end"), ":", convert=TRUE)
  LocMatrix <- LocMatrix %>% filter(chr==cc)
  index_chr <- LocMatrix %>% pull(name)
  map <- map[index_chr, index_chr]
  
  ### Normalization
  Normalization = as.character(opt["normalize"])
  if(Normalization == "average"){
    d <- nrow(map)
    map <- map / sum(map, na.rm=TRUE) * d * d
  }else if(Normalization == "probability"){
    map <- map / sum(map, na.rm=TRUE)
  }
  
  
  ### distance normalization
  MAP_RESOLUTION <- LocMatrix %>% head(n=1) %>% mutate(reso=end-start+1) %>% pull(reso)
  NUM_LINE <- nrow(map)
  MAXIMUM_DISTANCE <- D_location %>% filter(chr == cc) %>% mutate(d=end - start + 1) %>% arrange(desc(d)) %>% head(n=1) %>% pull(d)
  MAX_distance <- min(NUM_LINE-1, as.integer(MAXIMUM_DISTANCE / MAP_RESOLUTION))
  map.norm <- map
  for(d in 0:MAX_distance){
    index1 <- 1:(NUM_LINE - d)
    index2 <- index1 + d
    index3 <- cbind(LocMatrix[index1, "name"], LocMatrix[index2, "name"])
    index4 <- cbind(LocMatrix[index2, "name"], LocMatrix[index1, "name"])
    Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
    Score <- as.numeric(map[index3]) / Average
    map.norm[index3] <- Score
    map.norm[index4] <- Score
    rm(index1, index2, index3, index4, Average, Score)
  }
  
  D_table <- D_location %>% filter(chr == cc)
  G_LocMatrix <- GRanges(LocMatrix)
  
  Output_Map <- function(i){
    ### baitの部分のbinsを抽出
    D_bait <- D_table[i, c("bait_chr", "bait_start", "bait_end")]
    colnames(D_bait) <- c("chr", "start", "end")
    G_bait <- GRanges(D_bait)
    ov <- findOverlaps(G_bait, G_LocMatrix)
    index_bait <- G_LocMatrix$name[subjectHits(ov)]

    
    if(length(index_bait) > 1){
      D_score <- data.frame(LocMatrix, score=apply(map[index_bait,], 2, mean, na.rm=TRUE), score.norm=apply(map.norm[index_bait,], 2, mean, na.rm=TRUE))
    }else{
      D_score <- data.frame(LocMatrix, score=map[index_bait,], score.norm=map.norm[index_bait,])
    }
    
    D_score.extract <- D_score %>% filter(start >= D_table[i,"start"], end <= D_table[i, "end"])
    D_score.extract <- D_score.extract %>% select(chr, start, end, score, score.norm)
    
    FILE_out <- paste0(DIR_out, D_table[i,"name"], ".txt")
    write.table(D_score.extract, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  dummy <- sapply(1:nrow(D_table), Output_Map)
}




