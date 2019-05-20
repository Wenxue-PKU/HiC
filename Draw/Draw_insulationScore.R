# insulation scoreを描画する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"), default="NA", help="output png file"),
  make_option(c("--location"), default="NA", help="output location information as file"),
  make_option(c("--chr"),help="chromosome name"),
  make_option(c("--distance"), default="5", help="distance of bin start searching"),
  make_option(c("--window"), default="13", help="window size"),
  make_option(c("--threshold"), default="-0.3", help="threshold to define border"),
  make_option(c("--ylim"), default="1", help="y range"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--width"), default="NA", help="width of picture"),
  make_option(c("--height"), default="NA", help="height of picture")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_matrix <- as.character(opt["in"])
FILE_object <- sub("matrix", "rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
MapResolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1


CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(as.character(opt["end"]) == "all"){
  SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
  END <- max(as.numeric(LocMatrix[SameChromosome,3]));
}else{
  END <- as.numeric(as.character(opt["end"]))
}
IntraChromosome <- which(as.character(LocMatrix[,1]) == CHR)
IntraLoc <- LocMatrix[IntraChromosome,]
Region <- which((as.numeric(IntraLoc[,2]) >= START) & (as.numeric(IntraLoc[,3]) <= END))

map <- map[IntraChromosome, IntraChromosome]

# #============================================
# # Distance normalize
# #============================================
# NUM_LINE <- nrow(map)
# map_expect <- map
# for(d in 2:(NUM_LINE-1)){
#   index1 <- 1:(NUM_LINE - d)
#   index2 <- index1 + d
#   index3 <- cbind(index1, index2)
#   index4 <- cbind(index2, index1)
#   Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
#   map_expect[index3] <- Average
#   map_expect[index4] <- Average
# }
# map <- map / map_expect

#============================================
# Insulation score
#============================================
W <- as.numeric(as.character(opt["window"]))
D <- as.numeric(as.character(opt["distance"]))
InsulationSore <- function(i){
  sum(map[(i+D):(i+D+W-1),(i-D-W+1):(i-D)], na.rm=TRUE)
}
IS <- sapply((W + D):(nrow(map) - W - D +1), InsulationSore)
IS <- c(rep(NA, W+D-1), IS, rep(NA, W+D-1))
IS <- log2(IS / mean(IS, na.rm = TRUE))

IS <- IS[Region]



# boundaryを計算
threshold <- as.numeric(as.character(opt["threshold"]))
boundaries = c()
for(i in (W):(length(IS)-W)){
  condition <- IS[i] < threshold && IS[i] <= min(IS[(i-W):(i+W)], na.rm = TRUE)
  if(!is.na(condition) & condition){
    boundaries <- c(boundaries, i)
  }
}


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
  yrange <- as.numeric(as.character(opt["ylim"]))
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  plot(x3,IS, type='l', xaxs = "i", yaxs="i", axes=F, xlab='', ylab='', xlim=c(START, END), lwd=2, ylim=c(-yrange, yrange), bty = "n")
  abline(v=x3[boundaries], col='red')
  dummy <- dev.off()
}



### 座標を出力する
FILE_location <- as.character(opt["location"])
if(FILE_location != "NA"){
  bd <- rep(0, length(IS))
  bd[boundaries] = 1
  TADid <- cumsum(bd)
  TADid[boundaries] <- TADid[boundaries] - 0.5
  CountTADsize <- table(TADid)
  
  # domainのサイズは2以上で上位30%のISが0以上である
  TADorNot <- function(i){
    size <- CountTADsize[names(CountTADsize)==TADid[i]]
    score <- sort(IS[TADid==TADid[i]], decreasing = TRUE)
    num_30p <- as.integer(length(score) * 0.3)
    top_30_ave <- mean(score[1:num_30p])
    if( size > 1 && top_30_ave > 0){
      1
    }else{
      0
    }
  }
  
  TAD <- sapply(1:length(TADid), TADorNot)
  
  Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,2]) >= START) & (as.numeric(LocMatrix[,3]) <= END))
  OUTPUT <- cbind(LocMatrix[Region,], DI.borderCalc, DI, bd, TADid, TAD)
  write.table(OUTPUT, file=FILE_location, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE, append=FALSE)
}




