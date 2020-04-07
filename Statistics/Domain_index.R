#!/usr/bin/Rscript
# Domain score index


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-r", "--ref"), default="NA", help="Reference domain border file. (Use for checking location)"),
  make_option(c("-i", "--in"), default="NA", help="Target sample domain border file"),
  make_option(c("-m", "--matrix"), default="NA", help="matrix file"),
  make_option(c("--format"), default="rds", help="input file format (default: rds)"),
  make_option(c("-o", "--out"), default="NA", help="output file prefix")
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
options(scipen=10)


FILE_ref <- as.character(opt["ref"])
FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])

readDomain <- function(file){
  df <- fread(file, header=FALSE)
  if(ncol(df) == 7){
    colnames(df) <- c("chr", "start", "end", "border_strength", "border", "domain_num", "domain")
  }else if(ncol(df) == 8){
    colnames(df) <- c("chr", "start", "end", "border_strength", "border_strength_norm", "border", "domain_num", "domain")
  }
  df <- df %>% mutate(id=paste(chr, start, end, sep=":"))
  rownames(df) <- df %>% pull(id)
  df
}

D_ref <- readDomain(FILE_ref)
D_target <- readDomain(FILE_in)

FILE_format <- as.character(opt["format"])
FILE_matrix <- as.character(opt["matrix"])
if(FILE_format == "rds"){
  map <- readRDS(FILE_matrix)
}else if(FILE_format == "matrix"){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  cat("Unknown input file format")
  q()
}

r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1


#=============================================
# Average border strength
#=============================================
border <- D_ref %>% filter(border == 1) %>% pull(id)
ave_border <- mean(D_target[border, "border_strength"])


#=============================================
# Intr & inter domain (1Mb 以下)
#=============================================
D_table <- NULL
NUM_LINE <- nrow(map)
for(d in 2:(as.integer(500000 / Resolution))){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  Average <- mean(as.numeric(map[index3]), na.rm=TRUE, trim=0.1)
  
  index1 <- r[index1]
  index2 <- r[index2]
  
  df <- data.frame(doN1=as.numeric(D_ref[index1,"domain_num"]), doN2=as.numeric(D_ref[index2,"domain_num"]), 
                   doB1=as.numeric(D_ref[index1,"border"]), doB2=as.numeric(D_ref[index2,"border"]),
                   score=as.numeric(map[index3])/Average, stringsAsFactors = FALSE)
  
  df <- df %>% filter(!is.na(doN2 - doN1)) %>% filter(doB1==0, doB2==0, doN2 - doN1 < 2)  %>% mutate(category=if_else(doN2 - doN1 == 0, "intra", "inter"))
  D_table <- rbind(D_table , df %>% select(category, score))
}
D_table <- D_table %>% filter(!is.na(score))

# D_table %>% group_by(category) %>% summarize(n=n(), mean(score))


ave_intra <- D_table %>% filter(category=="intra") %>% pull(score)
ave_inter <- D_table %>% filter(category=="inter") %>% pull(score)

D_out <- data.frame(border_strength=ave_border, ratio=mean(ave_intra, na.rm=TRUE, trim=0.01) / mean(ave_inter, na.rm=TRUE, trim=0.01))
write.table(D_out, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




