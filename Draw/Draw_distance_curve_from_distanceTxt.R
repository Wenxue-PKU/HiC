# distance curveをmatrix からではなく、distanceファイルからデータを読み取って計算する


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--data"),help="Data directory"),
  make_option(c("-n", "--names"),help="name of target(s), separated by ,"),
  make_option(c("-o", "--out"),default="NA", help="output png file"),
  make_option(c("-c", "--chr"),default="chr1", help="target chromosome to calculate (default=chr1)"),
  make_option(c("--ratio"), default="FALSE", help="plotting log2 ratio from first data"),
  make_option(c("--width"), default="NA", help="width of output figure"),
  make_option(c("--height"), default="NA", help="height of output figure"),
  make_option(c("--xmin"), default="NA", help="minimum distance"),
  make_option(c("--xmax"), default="NA", help="maximum distance"),
  make_option(c("--ymin"), default="NA", help="minimum score"),
  make_option(c("--ymax"), default="NA", help="maximum score"),
  make_option(c("--color"), default="NA", help="color lists if not default"),
  make_option(c("--color2"), default="NA", help="fitting line colors"),
  make_option(c("--fill"), default="FALSE", help="TRUE if not axis"),
  make_option(c("-t", "--title"), default="", help="title ")
)
opt <- parse_args(OptionParser(option_list=option_list))

DIR <- as.character(opt["data"])
if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
  DIR <- paste(DIR, "/", sep="")
}

SAMPLES <- unlist(strsplit(as.character(opt["names"]), ","))
SAMPLE_NUMBER <- length(SAMPLES)
CHR <- as.character(opt["chr"])
TITLE <- as.character(opt["title"])

if(as.character(opt["color"]) == "NA"){
  c <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
}else{
  c <- unlist(strsplit(as.character(opt["color"]), ","))
}
c2 <- c
c <- adjustcolor(c, alpha.f = 0.2) 


if(as.character(opt["color2"]) != "NA"){
  c2 <- unlist(strsplit(as.character(opt["color2"]), ","))
}


DATA <- list()
for(NAME in SAMPLES){
  file <- paste(DIR, NAME, "_distance.txt", sep="")
  D <- read.table(file, header=T, sep="\t", stringsAsFactors = F)
  index <- which(D[,1] == CHR & D[,2] == CHR)
  DATA[[NAME]] <- D[index, c("distance", "probability")]
  colnames(DATA[[NAME]]) <- c("x", "y")
}

commonDis <- as.numeric(DATA[[1]][,"x"])
if(SAMPLE_NUMBER > 1){
  for(i in 2:SAMPLE_NUMBER){
    commonDis <- intersect(commonDis, as.numeric(DATA[[i]][,"x"]))
  }
}

for(NAME in SAMPLES){
  index <- DATA[[NAME]][,"x"] %in% commonDis
  # cat(NAME, sum(index), "\n")
  DATA[[NAME]] <- DATA[[NAME]][index, ]
}


# 1つめのファイルとのlog2 ratioを計算する
if(eval(parse(text=as.character(opt["ratio"])))){
  DATA.tmp <- DATA
  for(i in 2:SAMPLE_NUMBER){
    DATA.tmp[[i-1]][,"y"] <- log2(as.numeric(DATA[[i]][,"y"])/as.numeric(DATA[[1]][,"y"]))
  }
  SAMPLE_NUMBER <- SAMPLE_NUMBER -1
  DATA <- DATA.tmp
  rm(DATA.tmp)
}

if(as.character(opt["height"])=="NA"){
  h <- 20
}else{
  h <- as.numeric(as.character(opt["height"]))
}
if(as.character(opt["width"])=="NA"){
  w <- 20
}else{
  w <- as.numeric(as.character(opt["width"]))
}


FILE_OUT <- as.character(opt["out"])



if(as.character(opt["xmin"]) == "NA"){
  xmin <- min(commonDis, na.rm = TRUE)
}else{
  xmin <- as.numeric(as.character(opt["xmin"]))
}

if(as.character(opt["xmax"]) == "NA"){
  xmax <- max(commonDis, na.rm = TRUE)
}else{
  xmax <- as.numeric(as.character(opt["xmax"]))
}

if(as.character(opt["ymin"]) == "NA"){
  ymin <- min(as.numeric(DATA[[1]][["y"]]), na.rm = TRUE) * 0.8
}else{
  ymin <- as.numeric(as.character(opt["ymin"]))
}

if(as.character(opt["ymax"]) == "NA"){
  ymax <- max(as.numeric(DATA[[1]][["y"]]), na.rm = TRUE) * 1.2
}else{
  ymax <- as.numeric(as.character(opt["ymax"]))
}


if(sum((grep("\\.eps$", FILE_OUT))) == 1){
  if(w > 50){
    w<- w / 2.54
    h <- h / 2.54
  }
  cex.axis=0.7
  cex.lab=0.8
  postscript(file=FILE_OUT, horizontal=FALSE, onefile=FALSE, paper="special", height=h, width=w, family="Helvetica")
}
if(sum((grep("\\.png$", FILE_OUT))) == 1){
  cex.axis=1.4
  cex.lab=1.5
  png(file=FILE_OUT, width=w, height=h, bg="white", res = 200, units = "cm")
}
if(sum((grep("\\.pdf$", FILE_OUT))) == 1){
  if(w > 10){
    w<- w / 2.54
    h <- h / 2.54
  }
  cex.axis=1.4
  cex.lab=1.5
  pdf(file = FILE_OUT, height=h, width=w, useDingbats=FALSE)
}


LogAxis <- "xy"
Ylabels <- "Probability"
if(eval(parse(text=as.character(opt["ratio"])))){
  Ylabels <- "Log2 ratio of probability"
  LogAxis <- "x"
}



makeSmoothLines <- function(x,y){
  nn <- lowess(log10(x), log10(y), f=0.05)
  data.frame(x=10**nn$x, y=10**nn$y)
}
makeSmoothLines2 <- function(x,y){
  nn <- lowess(log10(x), y, f=0.05)
  data.frame(x=10**nn$x, y=nn$y)
}

if(eval(parse(text=as.character(opt["fill"])))){
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  plot(commonDis, DATA[[1]][["y"]], type='l', xaxs = "i", yaxs="i", xlab="Distance", ylab=Ylabels, 
       col=c[1], log = LogAxis, cex.axis=cex.axis, cex.lab=cex.lab, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), bty = "n", main=TITLE)
  if(SAMPLE_NUMBER > 1){
    for(n in 2:SAMPLE_NUMBER){
      par(new=T)
      plot(commonDis, DATA[[n]][["y"]], type='l', xaxs = "i", yaxs="i", axes=F, xlab="", ylab="", col=c[n], 
           log = LogAxis, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), bty = "n")
    }
  }
  pp <- makeSmoothLines(commonDis, DATA[[1]][["y"]])
  par(new=T)
  plot(pp$x, pp$y, type='l', xaxs = "i", yaxs="i", xlab="", ylab="", 
       col=c2[1], log = LogAxis, cex.axis=cex.axis, cex.lab=cex.lab, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), bty = "n")
  if(SAMPLE_NUMBER > 1){
    for(n in 2:SAMPLE_NUMBER){
      par(new=T)
      pp <- makeSmoothLines(commonDis, DATA[[n]][["y"]])
      plot(pp$x, pp$y, type='l', xaxs = "i", yaxs="i", axes=F, xlab="", ylab="", col=c2[n], 
           log = LogAxis, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), bty = "n")
    }
  }
}else{
  par(oma=c(0,0,0,0), mar=c(5,5,1,1))
  plot(commonDis, DATA[[1]][["y"]], type='l', xaxs = "i", yaxs="i", xlab="Distance", ylab=Ylabels,
       col=c[1], log = LogAxis, cex.axis=cex.axis, cex.lab=cex.lab, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), main=TITLE)
  if(SAMPLE_NUMBER > 1){
    for(n in 2:SAMPLE_NUMBER){
      par(new=T)
      plot(commonDis, DATA[[n]][["y"]], type='l', xaxs = "i", yaxs="i", axes=F, xlab="", ylab="", col=c[n], 
           log = LogAxis, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    }
  }
  if(eval(parse(text=as.character(opt["ratio"])))){
    pp <- makeSmoothLines2(commonDis, DATA[[1]][["y"]])
  }else{
    pp <- makeSmoothLines(commonDis, DATA[[1]][["y"]])
  }
  par(new=T)
  plot(pp$x, pp$y, type='l', xaxs = "i", yaxs="i", xlab="", ylab="",
       col=c2[1], log = LogAxis, cex.axis=cex.axis, cex.lab=cex.lab, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax))
  if(SAMPLE_NUMBER > 1){
    for(n in 2:SAMPLE_NUMBER){
      par(new=T)
      if(eval(parse(text=as.character(opt["ratio"])))){
        pp <- makeSmoothLines2(commonDis, DATA[[n]][["y"]])
      }else{
        pp <- makeSmoothLines(commonDis, DATA[[n]][["y"]])
      }
      plot(pp$x, pp$y, type='l', xaxs = "i", yaxs="i", axes=F, xlab="", ylab="", col=c2[n], 
           log = LogAxis, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    }
  }
}

LEGEND_label <- SAMPLES
if(eval(parse(text=as.character(opt["ratio"])))){
  LEGEND_label <- SAMPLES[-1]
}
legend("bottomleft", legend=LEGEND_label, lwd=2, lty=1, col=c2[1:SAMPLE_NUMBER], bty='n')

if(eval(parse(text=as.character(opt["ratio"])))){
  abline(h=0, col='blue', lty=3, lwd=2)
}
gabage <- dev.off()






