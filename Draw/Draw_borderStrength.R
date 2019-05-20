# directionary indexを描画する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"), default="NA", help="output png file"),
  make_option(c("--location"), default="NA", help="output location information as file"),
  make_option(c("--chr"),help="chromosome name"),
  make_option(c("--window"), default="13", help="window size"),
  make_option(c("--threshold"), default="0", help="window size"),
  make_option(c("--min_close"), default="2", help="window size"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--refresh"), default="FALSE", help="FALSE: use existing map object if eixts"),
  make_option(c("--width"), default="NA", help="width of picture"),
  make_option(c("--height"), default="NA", help="height of picture")
)
opt <- parse_args(OptionParser(option_list=option_list))



windowSize = as.numeric(as.character(opt["window"]))
BorderStrength <- function(m){
  minIndex <- windowSize + 1
  maxIndex <- dim(m)[1] - windowSize
  DI <- c()
  for(i in minIndex:maxIndex){
    # intra domain 前半
    A = 0
    for(j in (i-windowSize):(i-1)){
      for(p in (j+1):i){
        if(!is.na(m[j,p])){
          A = A + m[j,p]
        }
      }
    }
    # intra domain 後半
    B = 0
    for(j in i:(i+windowSize-1)){
      for(p in (j+1):(i+windowSize)){
        if(!is.na(m[j,p])){
          B = B + m[j,p]
        }
      }
    }
    # inter domain
    C = sum(m[(i-windowSize):(i-1),(i+1):(i+windowSize)], na.rm=TRUE)
    if(C==0){
      DI <- c(DI, 1)
    }else{
      DI <- c(DI, ((A+B)/C))
    }
  }
  DI
}

Flag_refresh <- eval(parse(text=as.character(opt["refresh"])))

FILE_matrix <- as.character(opt["in"])
FILE_object <- sub("matrix", "rds", FILE_matrix)
if(!file.exists(FILE_object) || Flag_refresh){
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

# DIを計算
DI <- BorderStrength(map[IntraChromosome , IntraChromosome])
DI.borderCalc <- c(rep(-99, windowSize), DI, rep(-99, windowSize))


# 大きすぎる値を整形する
maximumValue <- mean(DI) + sd(DI)*2
DI <- ifelse(DI > maximumValue, maximumValue, DI)
averageValue <- mean(DI)
DI <- c(rep(averageValue, windowSize), DI, rep(averageValue, windowSize))

DI <- as.numeric(scale(DI, scale = FALSE))
DI <- DI[Region]
DI.borderCalc <- DI.borderCalc[Region]
Xlength <- length(DI)


# boundaryを計算
check_bin <- as.numeric(as.character(opt["min_close"]))
threshold <- as.numeric(as.character(opt["threshold"]))
boundaries = c()
for(i in (check_bin+1):(Xlength-check_bin)){
  if(DI[i] > threshold && DI.borderCalc[i] >= max(DI.borderCalc[(i-check_bin):(i+check_bin)])){
    boundaries <- c(boundaries, i)
  }
}

FILE_OUT <- as.character(opt["out"])
if(FILE_OUT != "NA"){

  if(as.character(opt["height"])=="NA"){
    h <- 50
  }else{
    h <- as.numeric(as.character(opt["height"]))
  }
  if(as.character(opt["width"])=="NA"){
    w <- 1000
  }else{
    w <- as.numeric(as.character(opt["width"]))
  }
  
  
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    if(w > 100){
      w<- w / 72
      h <- h / 72
    }
    postscript(file=FILE_OUT, horizontal=FALSE, onefile=FALSE, paper="special", height=h, width=w, family="Helvetica")
  }
  if(sum((grep("\\.png$", FILE_OUT))) == 1){
    png(filename=FILE_OUT, width=w, height=h)
  }
  if(sum((grep("\\.pdf$", FILE_OUT))) == 1){
    if(w > 100){
      w<- w / 72
      h <- h / 72
    }
    pdf(file = FILE_OUT, height=h, width=w, useDingbats=FALSE)
  }
  
  x1 <- as.numeric(LocMatrix[Region,2])
  x2 <- as.numeric(LocMatrix[Region,3])
  x3 <- (x1 + x2) /2
  drawPoly <- function(i){
    polygon(c(x1[i], x2[i], x2[i], x1[i]), c(DI[i], DI[i], 0,0), col=ifelse(DI[i] > 0, 'orange', 'blue'), 
            border=NA, ylim=c(-1,3), xlim=c(START, END))
  }
  
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  plot(x1,DI, type='n', xaxs = "i", yaxs="i", axes=F, xlab='', ylab='', ylim=c(-1,3), xlim=c(START, END))
  dummy <- sapply(1:length(DI), drawPoly)  
  abline(v=x3[boundaries], col='red')
  grabage <- dev.off()
}

### 座標を出力する
FILE_location <- as.character(opt["location"])
if(FILE_location != "NA"){
  bd <- rep(0, length(DI))
  bd[boundaries] = 1
  TADid <- cumsum(bd)
  TADid[boundaries] <- TADid[boundaries] - 0.5
  CountTADsize <- table(TADid)
  
  # domainのサイズは2以上
  TADorNot <- function(i){
    size <- CountTADsize[names(CountTADsize)==TADid[i]]
    score <- sort(DI[TADid==TADid[i]])
    num_20p <- as.integer(length(score) * 0.3)
    bottom_20_ave <- mean(score[1:num_20p])
    if( size > 1 && bottom_20_ave < 0){
      1
    }else{
      0
    }
  }
  
  TAD <- sapply(1:length(TADid), TADorNot)
  
  Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,2]) >= START) & (as.numeric(LocMatrix[,3]) <= END))
  OUTPUT <- cbind(LocMatrix[Region,], format(DI.borderCalc, digits = 4), format(DI, digits = 4), bd, TADid, TAD)
  write.table(OUTPUT, file=FILE_location, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE, append=FALSE)
  
}





