#!/usr/bin/Rscript
# diff HiC

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="directories of map rds separated by ,"),
  make_option(c("-o", "--out"), default="NA", help="output file of p-values"),
  make_option(c("-g", "--group"), default="NA", help="groups (1 or 2) separated by ,. logFC >0 means 1<2")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(diffHic)))
suppressWarnings(suppressMessages(library(edgeR)))
FILE_samples <- as.character(opt["in"])
FILE_samples <- unlist(strsplit(FILE_samples, ","))
group <- as.character(opt["group"])
group <- factor(unlist(strsplit(group, ",")))
FILE_output <- as.character(opt["out"])

if(file.exists(FILE_output)){
  file.remove(FILE_output)
}

SAMPLES <- 1:length(FILE_samples)

for(chr in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")){
  
  # cat("start ", chr, "...\n")
  
  getMap <- function(dir){
    if(substring(dir, nchar(dir), nchar(dir)) != "/"){
      dir <- paste(dir, "/", sep="")
    }
    file <- paste0(dir, "/Raw/", chr, ".rds")
    readRDS(file)
  }
  
  # cat("Load maps ...\n")
  r.common <- c()
  maps <- list()
  for(ss in SAMPLES){
    maps[[ss]] <- getMap(FILE_samples[ss])
    if(length(r.common) == 0){
      r.common <- rownames(maps[[ss]])
    }else{
      r.common <- intersect(r.common, rownames(maps[[ss]]))
    }
  }
  LocList <- strsplit(r.common, ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  df <- data.frame(chr=as.character(LocMatrix[,1]), start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]))
  regions <- makeGRangesFromDataFrame(df)
  rm(LocList, df)
  
  
  # cat("Convert map to contactMatrix ...\n")
  lib.size <- c()
  cm <- list()
  for(ss in SAMPLES){
    map.filter <- maps[[ss]][r.common, r.common]
    cm[[ss]] <- ContactMatrix(map.filter, regions, regions)
    lib.size <- c(lib.size, sum(map.filter[upper.tri(map.filter, diag = TRUE)], na.rm = TRUE))
  }
  rm(maps)
  
  to.keep <- as.matrix(cm[[1]])
  for(i in 2:length(cm)){
    to.keep <- to.keep + as.matrix(cm[[i]])
  }
  to.keep <- to.keep!=0
  
  
  for(i in 1:length(cm)){
    iset <- deflate(cm[[i]], extract=to.keep)
    if(i==1){
      data <- iset
    }else{
      data <- cbind(data, iset)
    }
  }
  
  # cat("Make data ...\n")
  interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
  data$totals <- lib.size
  
  # cat("edgeR analysis ...\n")
  keep <- aveLogCPM(asDGEList(data, group=group)) > 0
  data <- data[keep,]
  data <- normOffsets(data, type="loess", se.out=TRUE)
  y <- asDGEList(data)
  
  # design matrixの作成
  design <- model.matrix(~group)
  
  
  # cat("Estimate dispersion...\n")
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)
  
  
  result <- glmQLFTest(fit)
  adj.p <- p.adjust(result$table$PValue, method="BH")
  NUM <- sum(result$table$PValue <= 0.05)
  
  useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
  inter.frame <- as.data.frame(interactions(data))[,useful.cols]
  results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
  results.r <- results.r[order(results.r$PValue),]
  if(NUM == 0){
    results.r <- results.r[1:10,]
  }else{
    results.r <- results.r[results.r$PValue < 0.05,]
  }
  if(chr == "chr1"){
    write.table(results.r, file=FILE_output, sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, append = TRUE)
  }else{
    write.table(results.r, file=FILE_output, sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE, append = TRUE)
  }
}


