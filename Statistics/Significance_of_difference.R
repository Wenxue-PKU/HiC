#!/usr/bin/Rscript
# Define significant pairs by wilcox test

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="all data table"),
  make_option(c("--data"), default="NA", help="Data directory"),
  make_option(c("--chr"), default="NA", help="chromosome"),
  make_option(c("--con1"), default="NA", help="first sample groups. Separated by ,"),
  make_option(c("--con2"), default="NA", help="second sample groups. Separated by ,"),
  make_option(c("-o", "--out"), default="NA", help="output table")
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(stringr)))

CHR <- as.character(opt["chromosome"])

FILE_data <- as.character(opt["in"])
DATA_head <- read.table(FILE_data, header=TRUE, nrows = 5, stringsAsFactors = FALSE)
classes <- sapply(DATA_head, class)
D_data <- read.table(FILE_data, header=TRUE, colClasses = classes, stringsAsFactors = FALSE)
rm(DATA_head, classes)

SAMPLEs <- D_data %>% select(-Loc1, -Loc2) %>% colnames()


Left <- str_split_fixed(D_data$Loc1, ":", 3)
Right <- str_split_fixed(D_data$Loc2, ":", 3)
Distance <- (as.numeric(Right[,2])+as.numeric(Right[,3]))/2 - (as.numeric(Left[,2])+as.numeric(Left[,3]))/2
D_data <- D_data %>% mutate(distance=Distance)
rm(Left, Right, Distance, FILE_data)


DIR_data <- as.character(opt["data"])
if(substring(DIR_data, nchar(DIR_data), nchar(DIR_data)) != "/"){
  DIR_data <- paste0(DIR_data, "/")
}

getDistance <- function(name){
  FILE_dis <- paste0(DIR_data, name, "_distance.txt")
  DATA_head <- read.table(FILE_dis, header=TRUE, nrows = 5, stringsAsFactors = FALSE)
  classes <- sapply(DATA_head, class)
  D_dis <- read.table(FILE_dis, header=TRUE, colClasses = classes, stringsAsFactors = FALSE)
  rm(DATA_head)
  name_con <- paste0(name, "_background")
  D_dis <- D_dis %>% filter(chromosome1 == CHR & chromosome2 == CHR) %>% select(distance, average) %>% rename(!!as.name(name_con) := average)
  setDT(D_data)
  setDT(D_dis)
  setkey(D_data, distance)
  setkey(D_dis, distance)
  D_data <<- D_dis[D_data, roll="nearest"]
}

for(ss in SAMPLEs){
  getDistance(ss)
}


for(name in SAMPLEs){
  name_fd <- paste0(name, "_fc")
  name_con <- paste0(name, "_background")
  D_data <- D_data %>% mutate(!!as.name(name_fd) := !!as.name(name) / !!as.name(name_con))
}

### P-value
MultiPval <- function(c1, c2){
  getPval <- function(i){
    aa <- D_data[i,c1]
    bb <- D_data[i,c2]
    p <- wilcox.test(aa, bb, alternative = "two.sided")  
    p$p.value
  }
  p <- sapply(1:nrow(D_data), getPval)
  p
}

group1 <- unlist(strsplit(as.character(opt["con1"]), ","))
group2 <- unlist(strsplit(as.character(opt["con2"]), ","))


D_data <- D_data %>% mutate(Pval=MultiPval(group1, group2)) %>% mutate(FDR=p.adjust(Pval, method = "BH"))

FILE_out <- as.character(opt["out"])
write.table(FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)







