#!/usr/bin/Rscript
# Big domain score index


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="matrix file"),
  make_option(c("--format"), default="rds", help="input file format (default: rds)"),
  make_option(c("-o", "--out"), default="NA", help="output file"),
  make_option(c("-s", "--size"), default="big", help="'big' domain or 'small' domains"),
  make_option(c("--location"), default="FALSE", help="Output location information or not (default:FALSE)")
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
options(scipen=10)


DOMAIN_SIZE <- as.character(opt["size"])
FLAG_location <- eval(parse(text=as.character(opt["location"])))

if(DOMAIN_SIZE == "big"){
  FILE_ref <- "G:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-04-12_domain_statistics/time3_r1/Manual_defined_domain.txt"
  FILE_ref <- "~/Mount/Genomics/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-04-12_domain_statistics/time3_r1/Manual_defined_domain.txt"
}else{
  FILE_ref <- "G:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-04-12_domain_statistics/WT1/Small_ALL.txt"
  FILE_ref <- "~/Mount/Genomics/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-04-12_domain_statistics/WT1/Small_ALL.txt"
}


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

FILE_format <- as.character(opt["format"])
FILE_matrix <- as.character(opt["in"])
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
# Intr & inter domain 
#=============================================
#---------------------------------------------
# Threhold
#---------------------------------------------
if(DOMAIN_SIZE == "big"){
  ### Bigdomain (100kb 以上 500kb 以下)
  T_MIN <- 100000
  T_MAX <- 500000
}else if(DOMAIN_SIZE == "small"){
  ### Smalldomain (10kb 以上 150kb 以下)
  T_MIN <- 10000
  T_MAX <- 150000
}else{
  cat("Please select big or small domain")
  q()
}

D_table <- NULL
NUM_LINE <- nrow(map)
for(d in (as.integer(T_MIN / Resolution)):(as.integer(T_MAX / Resolution))){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  Average <- mean(as.numeric(map[index3]), na.rm=TRUE, trim=0.1)
  
  index1 <- r[index1]
  index2 <- r[index2]
  
  if(FLAG_location){
    df <- data.frame(doNL=as.numeric(D_ref[index1,"domain_num"]), doNR=as.numeric(D_ref[index2,"domain_num"]),
                     Raw=as.numeric(map[index3]),
                     NormScore=as.numeric(map[index3])/Average, 
                     id1=index1, id2=index2,
                     stringsAsFactors = FALSE)
  }else{
    df <- data.frame(doNL=as.numeric(D_ref[index1,"domain_num"]), doNR=as.numeric(D_ref[index2,"domain_num"]),
                     Raw=as.numeric(map[index3]),
                     NormScore=as.numeric(map[index3])/Average, stringsAsFactors = FALSE)
  }

  
  df <- df %>% filter(!is.na(NormScore)) %>% filter(!is.na(doNL), !is.na(doNR))
  D_table <- rbind(D_table , df)
}

write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




