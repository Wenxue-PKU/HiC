### ChIA-PETのような相互作用のグラフを作成する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"), help="outpu file name png or eps"),
  make_option(c("--color"), default="matlab", help="color matlab or gentle"),
  make_option(c("--chr"),help="chromosome name"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--width"), default="500", help="width of picture"),
  make_option(c("--height"), default="250", help="height of picture"),
  make_option(c("--mininteraction"), default="NA", help="minimum interaction to draw"),
  make_option(c("--topxxx"), default="NA", help="draw only top xxx lines"),
  make_option(c("--mindistance"), default="40000", help="minimum distance (bp) to draw")
)
opt <- parse_args(OptionParser(option_list=option_list))

### 3Dspline.Rを読み込む
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
other.name <- paste(sep="/", script.basename, "3Dspline.R")
source(other.name)

c <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                        "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))( 101 )
if(as.character(opt["color"]) == "gentle"){
  suppressWarnings(suppressMessages(library("RColorBrewer")))
  c <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(101))  # gentle color
}
if(as.character(opt["color"]) == "red"){
  suppressWarnings(library("RColorBrewer"))
  red <- brewer.pal(9, "Reds")[-1]
  c <- colorRampPalette(red)(101)  # gentle color
}


map <- as.matrix(read.table(as.character(opt["in"])), header=TRUE, check.names = FALSE)
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(as.character(opt["end"]) == "all"){
  SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
  END <- max(as.numeric(LocMatrix[SameChromosome,3]));
}else{
  END <- as.numeric(as.character(opt["end"]))
}
Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,3]) >= START) & (as.numeric(LocMatrix[,2]) <= END))
map <- map[Region, Region]
LocMatrix <- LocMatrix[Region,]


# matrixの解像度
resolution <- as.numeric(LocMatrix[2,2]) - as.numeric(LocMatrix[1,2])


# 値でfilterをかける
if(as.character(opt["mininteraction"]) != "NA"){
  map <- ifelse(map < as.numeric(as.character(opt["mininteraction"])), NA, map)
}

# ２点の間の距離でfilterをかける
minDistance <- as.numeric(as.character(opt["mindistance"]))
minBins <- as.integer(minDistance / resolution)
filtMinDistance <-function(i,j){
  map[i,i+j] <<- NA
}
for(i in 0:minBins){
  dummy <- mapply(filtMinDistance, 1:(nrow(map)-i), i)
}


### 下三角をすべて除去する（上三角と同じなので）
map[lower.tri(map)] <- NA


FILE_OUT <- as.character(opt["out"])
if(sum((grep("\\.eps$", FILE_OUT))) == 1){
  postscript(file=FILE_OUT, horizontal=FALSE, onefile=FALSE, paper="special", height =2.2, width=4)
}
if(sum((grep("\\.png$", FILE_OUT))) == 1){
  png(filename=FILE_OUT, width=as.numeric(as.character(opt["width"])), height=as.numeric(as.character(opt["height"])), units="px", bg="white")
}

par(oma=c(0,0,0,0), mar=c(0,0,0,0))
Graph_heigth <- (END - START)/2
plot(c(START, END), c(0, Graph_heigth), type='n', 
     xlab="", ylab="", xaxs="i", yaxs="i", axes=FALSE, bty="n", ylim=c( Graph_heigth, 0), xlim=c(START, END))

### 描画する組み合わせのindex
DrawIndex <- which(!is.na(map), arr.ind = TRUE)
DrawIndex <- cbind(DrawIndex, mapply(function(a, b) map[a,b], DrawIndex[,1], DrawIndex[,2]))
DrawIndex <- DrawIndex[order(DrawIndex[,3], decreasing = TRUE),]

# for(i in 1:10){
#   cat(paste(LocMatrix[DrawIndex[i, 1],2], LocMatrix[DrawIndex[i, 2],2], DrawIndex[i, 3]), "\n")
# }

#cat(nrow(DrawIndex))

DrawNumbers <- nrow(DrawIndex)
if(as.character(opt["topxxx"]) != "NA"){
  DrawNumbers <- as.numeric(as.character(opt["topxxx"]))
}

### 最大値と最小値を調べる
MaxScore <- max(DrawIndex[1:DrawNumbers, 3], na.rm = TRUE)
MinScore <- min(DrawIndex[1:DrawNumbers, 3], na.rm = TRUE)

drawConnecLine <- function(i){
  k1 <- as.numeric(DrawIndex[i,1])
  k2 <- as.numeric(DrawIndex[i,2])
  
  p1 <- as.numeric(LocMatrix[k1, 2])
  p2 <- as.numeric(LocMatrix[k2, 2])
  score <- map[k1,k2]
    
  middle <- (p1 + p2)/2 + resolution/2
  radius <- abs(p2 - p1)/2
  Px <- c()
  Py <- c()
  for(i in seq(0,180,10)){
    radian <- pi / 180 * i
    x <- radius * cos(radian) + middle
    y <- radius * sin(radian)
    Px <- c(Px, x)
    Py <- c(Py, y)
  }
  out <- Spl2.int(Px, Py, 0, 10)
  
  score_percentage <- (score - MinScore) / (MaxScore - MinScore) * 100
  color <- c[round(score_percentage) + 1]

  lines(out[,1], out[,2], lwd=1, col=color)
}

dummy <- sapply(DrawNumbers:1, drawConnecLine)
dummy <- dev.off()



