# 指定した領域同士のスコアを出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("--chr"),help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--chr2"), default="NULL", help="chromosome2 name. defula is same to chr"),
  make_option(c("--start2"), default="NULL", help="start2 position. default is same to start"),
  make_option(c("--end2"), default="NULL", help="end2 position. default is same to end"),
  make_option(c("--method"), default="max", help="max, min, mean, median")
)
opt <- parse_args(OptionParser(option_list=option_list))


map <- as.matrix(read.table(as.character(opt["in"])), header=TRUE, check.names = FALSE)
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)


CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(CHR=="all"){
  Region <- 1:nrow(map)
}else{
  if(as.character(opt["end"]) == "all"){
    SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
    END <- max(as.numeric(LocMatrix[SameChromosome,3]));
  }else{
    END <- as.numeric(as.character(opt["end"]))
  }
  Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,3]) >= START) & (as.numeric(LocMatrix[,2]) <= END))
}
CHR2 <- as.character(opt["chr2"])
START2 <- as.numeric(as.character(opt["start2"]))
if(CHR2=="all"){
  Region2 <- 1:nrow(map)
}else{
  if(as.character(opt["end2"]) == "all"){
    SameChromosome <- which(as.character(LocMatrix[,1]) == CHR2)
    END2 <- max(as.numeric(LocMatrix[SameChromosome,3]));
  }else{
    END2 <- as.numeric(as.character(opt["end2"]))
  }
  Region2 <- which((as.character(LocMatrix[,1]) == CHR2) & (as.numeric(LocMatrix[,3]) >= START2) & (as.numeric(LocMatrix[,2]) <= END2))
}
map <- as.numeric(map[Region, Region2])


METHOD=as.character(opt["method"])
if(METHOD == "max"){
  out <- max(map, na.rm=TRUE)
}else if(METHOD == "min"){
  out <- min(map, na.rm=TRUE)
}else if(METHOD == "median"){
  out <- median(map, na.rm = TRUE)
}else if(METHOD == "mean"){
  out <- mean(map, na.rm = TRUE)
}else{
  cat("unknown method: ", METHOD, "\n", sep="")
}
#cat(map, "\n", sep="")
cat(out)



