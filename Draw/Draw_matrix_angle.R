#!/usr/bin/Rscript
# HiC matrixの描画のためのテキストファイルを出力する
# 横に寝かしたデータに対応


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="directory of matrices located or rds file"),
  make_option(c("-o", "--out"),help="output directory"),
  make_option(c("--location"), default="NA", help="file describe target loci. text or excel file. chr, start, end, name"),
  make_option(c("--max_distance"), default="20000000", help="maximum distance for hic "),
  make_option(c("--sufix"), default="NA", help="sufix for output file name."),
  make_option(c("--normalize"), default="NA", help="NA, average: average will be 1, probability: score were divided by total read"),
  make_option(c("--na"), default="remove", help="how to treat na value. min, na, zero, remove. min replace with minimum value. zero replace to zero")
)
opt <- parse_args(OptionParser(option_list=option_list))

options(scipen=10)
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))


# ### test data
FILE_location <- "T:/Project/018_20191128_Noma_Senescent3/out/2020-08-31_example_map_hic_at_EP/data/draw_location_test.bed"
DIR_out <- "T:/Project/018_20191128_Noma_Senescent3/out/2020-08-31_example_map_hic_at_EP/data/"
DIR_in <- "T:/Project/018_20191128_Noma_Senescent3/data/IMR90_G_bmix/40kb/ICE2"
MAXIMUM_DISTANCE <- 6000000
NAME_SUFIX <- "_IMR90_G_bmix"







### location file
# chr, start, end, name
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
MAXIMUM_DISTANCE <- as.numeric(as.character(opt["max_distance"]))
TARGET <- as.character(opt["target"])


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
  
  ### Convert to angled HIC
  MAP_RESOLUTION <- LocMatrix %>% head(n=1) %>% mutate(reso=end-start+1) %>% pull(reso)
  
  NUM_LINE <- nrow(map)
  MAX_distance <- min(NUM_LINE-1, as.integer(MAXIMUM_DISTANCE / MAP_RESOLUTION))
  D_map <- NULL
  for(d in 0:MAX_distance){
    index1 <- 1:(NUM_LINE - d)
    index2 <- index1 + d
    index3 <- cbind(LocMatrix[index1, "name"], LocMatrix[index2, "name"])
    Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
    
    Middle <- as.integer((LocMatrix[index1, "start"] + LocMatrix[index2, "end"] + 1)/2)
    df1 <- data.frame(distance=LocMatrix[index2, "start"] - LocMatrix[index1, "start"], score=map[index3])
    df2 <- rbind(df1 %>% mutate(location=Middle - MAP_RESOLUTION/4), df1 %>% mutate(location=Middle + MAP_RESOLUTION/4))
    df2 <- df2 %>% arrange(location)
    D_map <- rbind(D_map, df2)
    rm(index1, index2, index3, Average, Middle, df1, df2)
  }
  D_disAve <- D_map %>% filter(!is.na(score)) %>% group_by(distance) %>% summarize(ave=mean(score), .groups = 'drop') %>% as.data.frame()
  
  D_map <- dplyr::left_join(D_map, D_disAve, by="distance")
  D_map <- D_map %>% mutate(score.norm=score/ave)
  D_map <- D_map %>% select(location, distance, score, score.norm)
  
  D_table <- D_location %>% filter(chr == cc)
  
  Output_Map <- function(i){
    D_map.extract <- D_map %>% filter(location >= D_table[i,"start"], location <= D_table[i, "end"])
    
    # replace Na
    if(as.character(opt["na"]) == "min"){
      D_map.extract <- D_map.extract %>% mutate(score = if_else(is.na(score), Min, score))
    }else if(as.character(opt["na"]) == "zero"){
      D_map.extract <- D_map.extract %>% mutate(score = if_else(is.na(score), 0, score))
    }else if(as.character(opt["na"]) == "remove"){
      D_map.extract <- D_map.extract %>% filter(!is.na(score))
    }
    
    D_map.extract <- D_map.extract %>% mutate(chr=cc)
    D_map.extract <- D_map.extract %>% select(chr, location, distance, score, score.norm)
    
    FILE_out <- paste0(DIR_out, D_table[i,"name"], ".txt")
    write.table(D_map.extract, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
  dummy <- sapply(1:nrow(D_table), Output_Map)
  
}






