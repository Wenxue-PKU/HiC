# BED fileで示す場所を出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-o", "--out"), help="outpu file name png or eps"),
  make_option(c("--chr"),help="chromosome name"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="1000", help="end position. all for end of the chromosome"),
  make_option(c("--width"), default="NA", help="width of picture"),
  make_option(c("--height"), default="NA", help="height of picture"),
  make_option(c("--position"), default="bottom", help="bottom or top axis position"),
  make_option(c("--cairo"), default="TRUE", help="use cairo for output")
)
opt <- parse_args(OptionParser(option_list=option_list))

CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
END <- as.numeric(as.character(opt["end"]))
FLAG_cairo <- eval(parse(text=as.character(opt["cairo"])))

FILE_OUT <- as.character(opt["out"])
if(as.character(opt["height"])=="NA"){
  h <- 120
}else{
  h <- as.numeric(as.character(opt["height"]))
}
if(as.character(opt["width"])=="NA"){
  w <- 1414
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
  
  if(FLAG_cairo){
    suppressWarnings(suppressMessages(library(Cairo)))
    CairoPNG(FILE_OUT, width=w, height=h)
  }else{
    png(filename=FILE_OUT, width=w, height=h)
  }
}
if(sum((grep("\\.pdf$", FILE_OUT))) == 1){
  if(w > 100){
    w<- w / 72
    h <- h / 72
  }
  pdf(file = FILE_OUT, height=h, width=w)
}

position <- as.character(opt["position"])
if(position == "top"){
  par(oma=c(0,0,0,0), mar=c(0,0,6.5,0), ps=18)
  plot(c(START, END), c(0,1), type='n', xaxs = "i", yaxs="i", xlab="", ylab="", ylim=c(0,1), axes=FALSE, xlim=c(START, END))
  axis(3, cex.axis = 1.4, lwd=2, tck=-2, line=0, mgp=c(0,2,0))
  mtext(paste("Position in Chromosome ", CHR, " (bp)", sep=""), side=3, line=4, cex=1.3)
}else{
  par(oma=c(0,0,0,0), mar=c(6.5,0,0,0), ps=18)
  plot(c(START, END), c(0,1), type='n', xaxs = "i", yaxs="i", xlab="", ylab="", ylim=c(0,1), axes=FALSE, xlim=c(START, END))
  axis(1, cex.axis = 1.4, lwd=2, tck=-0.8, line=0, mgp=c(0,3,0))
  mtext(paste("Position in Chromosome ", CHR, " (bp)", sep=""), side=1, line=5.2, cex=1.3)
}
dummy <- dev.off()


