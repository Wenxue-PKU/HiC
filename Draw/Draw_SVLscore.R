#!/usr/bin/Rscript
# SVL score (short range / long range)を計算する


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output png file"),
  make_option(c("--threshold"),default="2000000", help="threshold of short"),
  make_option(c("--max"),default="50000000", help="maximum distance of long range associations"),
  make_option(c("--score"), default="NA", help="output SVL score"),
  make_option(c("--width"), default="NA", help="width of picture"),
  make_option(c("--height"), default="NA", help="height of picture"),
  make_option(c("--chr"),help="chromosome name"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--ymax"), default="NA", help="Maximum score of Y-axis"),
  make_option(c("--ymin"), default="NA", help="Minimum score of Y-axis")
)
opt <- parse_args(OptionParser(option_list=option_list))



FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
rownames(LocMatrix) <- r

THRESHOLD <- as.numeric(as.character(opt["threshold"]))
MAXIMUM_DISTANCE <- as.numeric(as.character(opt["max"]))


CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(as.character(opt["end"]) == "all"){
  SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
  END <- max(as.numeric(LocMatrix[SameChromosome,3]));
}else{
  END <- as.numeric(as.character(opt["end"]))
}
Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,2]) >= START) & (as.numeric(LocMatrix[,3]) <= END))


getSVLscore <- function(i){
  intra <- as.character(LocMatrix[i,1]) == as.character(LocMatrix[,1])
  distance <- abs(as.numeric(LocMatrix[i,2]) - as.numeric(LocMatrix[,2]))
  short <- intra & distance < THRESHOLD & distance > 0
  long <- intra & distance > THRESHOLD & distance < MAXIMUM_DISTANCE
  score <- log2(sum(map[i, r[long]], na.rm = TRUE) / sum(map[i, r[short]], na.rm = TRUE))
  score
}

SVL <- sapply(r, getSVLscore)



FILE_OUT <- as.character(opt["out"])
if(FILE_OUT != "NA"){
  if(as.character(opt["height"])=="NA"){
    if(sum((grep("\\.eps$", FILE_OUT))) == 1){
      h <- 0.8
    }else{
      h <- 50
    }
  }else{
    h <- as.numeric(as.character(opt["height"]))
  }
  if(as.character(opt["width"])=="NA"){
    if(sum((grep("\\.eps$", FILE_OUT))) == 1){
      w <- 4
    }else{
      w <- 500
    }
  }else{
    w <- as.numeric(as.character(opt["width"]))
  }
  
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    postscript(file=FILE_OUT, horizontal=FALSE, onefile=FALSE, paper="special", height=h, width=w, family="Helvetica")
  }
  if(sum((grep("\\.png$", FILE_OUT))) == 1){
    png(filename=FILE_OUT, width=w, height=h)
  }
  
  x1 <- as.numeric(LocMatrix[Region,2])
  x2 <- as.numeric(LocMatrix[Region,3])
  x3 <- (x1 + x2) /2
  SVL_region <- SVL[Region]
  
  if(as.character(opt["ymax"]) == "NA"){
    ymax <- max(SVL_region, na.rm=TRUE)
  }else{
    ymax <- as.numeric(as.character(opt["ymax"]))
  }
  if(as.character(opt["ymin"]) == "NA"){
    ymin <- min(SVL_region, na.rm=TRUE)
  }else{
    ymin <- as.numeric(as.character(opt["ymin"]))
  }
  
  drawPoly <- function(i){
    polygon(c(x1[i], x2[i], x2[i], x1[i]), c(SVL_region[i], SVL_region[i], 0,0), col=ifelse(SVL_region[i] > 0, 'orange', 'blue'), 
            border=NA, xlim=c(START, END), ylim=c(ymin,ymax))
  }
  
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  plot(x1,SVL_region, type='n', xaxs = "i", yaxs="i", axes=F, xlab='', ylab='', ylim=c(ymin,ymax), xlim=c(START, END))
  dummy <- sapply(1:length(SVL_region), drawPoly)
  trash <- dev.off()
}

### スコアを出力する
FILE_score <- as.character(opt["score"])
if(FILE_score != "NA"){
  OUTPUT <- cbind(LocMatrix, SVL)
  write.table(OUTPUT, file=FILE_score, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE, append=FALSE)
  
}
