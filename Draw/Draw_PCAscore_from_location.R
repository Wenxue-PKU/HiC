#!/usr/bin/Rscript
# PCAの計算済みのlocationファイルを読み取ってグラフを作成する


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="location file from PCA"),
  make_option(c("-o", "--out"), default="NA", help="output png file"),
  make_option(c("--chr"),help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--ymax"), default="NA", help="Maximum score of Y-axis"),
  make_option(c("--ymin"), default="NA", help="Minimum score of Y-axis"),
  make_option(c("--width"), default="NA", help="width of output figure"),
  make_option(c("--line"), default="NA", help="line of plotting"),
  make_option(c("--fill"), default="FALSE", help="fill image or not"),
  make_option(c("--title"), default="", help="title of the graph"),
  make_option(c("--height"), default="NA", help="height of output figure")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_OUT <- as.character(opt["out"])
DATA <- read.table(as.character(opt["in"]), sep="\t", header=F, stringsAsFactors = FALSE, check.names = FALSE)


CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(as.character(opt["end"]) == "all"){
  END <- max(as.numeric(DATA[as.character(DATA[,1])==CHR,3]))
}else{
  END <- as.numeric(as.character(opt["end"]))
}


if(as.character(opt["height"])=="NA"){
  h <- 100
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
  png(file=FILE_OUT, width=w, height=h, units="px", bg="white")
}
if(sum((grep("\\.pdf$", FILE_OUT))) == 1){
  if(w > 100){
    w<- w / 72
    h <- h / 72
  }
  pdf(file = FILE_OUT, height=h, width=w, useDingbats=FALSE)
}

NUM <- sort(as.numeric(DATA[,4]), na.last = NA)
if(as.character(opt["ymax"]) == "NA"){
  YMAX <- NUM[round(length(NUM)*0.95+1)]
}else{
  YMAX=as.numeric(as.character(opt["ymax"]))
}
if(as.character(opt["ymin"]) == "NA"){
  YMIN <- NUM[length(NUM)*0.05+1]
}else{
  YMIN=as.numeric(as.character(opt["ymin"]))
}

index <- which(as.character(DATA[,1]) == CHR & as.numeric(DATA[,3]) >= START &  as.numeric(DATA[,2]) <= END)
x1 <- as.numeric(DATA[index,2])
x2 <- as.numeric(DATA[index,3])
y <- as.numeric(DATA[index,4])
compartment <- as.character(DATA[index,ncol(DATA)])

drawPoly <- function(i){
  # polygon(c(x1[i],x2[i],x2[i],x1[i]), c(y[i], y[i], 0, 0), col=ifelse(y[i]>0,'red', 'blue'), ylim=c(YMIN,YMAX), xlim=c(START,END), border=NA)
  polygon(c(x1[i],x2[i],x2[i],x1[i]), c(y[i], y[i], 0, 0), col=ifelse(compartment[i]=="A",'red', 'blue'), ylim=c(YMIN,YMAX), xlim=c(START,END), border=NA)
}

if(eval(parse(text=as.character(opt["fill"])))){
  par(oma=c(0,0,0,0), mar=c(1,0,1,0))
  plot(c(START, END), c(YMIN,YMAX), type='n', xaxs = "i", yaxs="i", xlab="", ylab="", ylim=c(YMIN,YMAX), xlim=c(START,END), axes=FALSE, bty = "n")
  dummy <- sapply(1:length(y), drawPoly)
}else{
  par(oma=c(0,0,0,0), mar=c(3,4,1,1))
  plot(c(START, END), c(YMIN,YMAX), type='n', xaxs = "i", yaxs="i", xlab="", ylab="", ylim=c(YMIN,YMAX), xlim=c(START,END))
  dummy <- sapply(1:length(y), drawPoly)
}

if(as.character(opt["line"]) != "NA"){
  abline(h=as.numeric(as.character(opt["line"])), lty=3)
}
if(as.character(opt["title"]) != ""){
  legend("topright", legend=as.character(opt["title"]), bty='n', plot=T, seg.len = 0.01)
}


gabage <- dev.off()