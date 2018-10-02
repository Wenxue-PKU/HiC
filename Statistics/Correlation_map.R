# ２つのmapの相関を求める

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--file1"),help="matrix file1"),
  make_option(c("-b", "--file2"),help="matrix file2"),
  make_option(c("--out"), default="NA", help="output graph file"),
  make_option(c("--NoPlot"), default="FALSE", help="do not plot (draw only axis)"),
  make_option(c("--title1"), default="", help="title of sample"),
  make_option(c("--title2"), default="", help="title of sample"),
  make_option(c("--out_intra"), default="TRUE", help="output is intra-chromosome or not"),
  make_option(c("--min_dis"), default="0", help="minimum distance"),
  make_option(c("--max_dis"), default="9999999", help="maximum distance"),
  make_option(c("--fill"), default="FALSE", help="draw with no space"),
  make_option(c("--distance"), default="all", help="maximum distance considered for intra"),
  make_option(c("--width"), default="NA", help="width of picture"),
  make_option(c("--height"), default="NA", help="height of picture")
)
opt <- parse_args(OptionParser(option_list=option_list))

# options(scipen=10)
map1 <- as.matrix(read.table(as.character(opt["file1"]), header=TRUE, check.names = FALSE))
map2 <- as.matrix(read.table(as.character(opt["file2"]), header=TRUE, check.names = FALSE))
r1 <- rownames(map1)
r2 <- rownames(map2)
r.common <- intersect(r1, r2)
map1 <- map1[r.common, r.common]
map2 <- map2[r.common, r.common]

LocList <- strsplit(r.common, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

Chromosomes <- unique(LocMatrix[,1])

t.all <- cor(as.numeric(map1), as.numeric(map2), use = "pairwise.complete.obs")


MinDis <- as.numeric(as.character(opt["min_dis"]))
MaxDis <- as.numeric(as.character(opt["max_dis"]))

map1.intra <- c()
map2.intra <- c()
for(i in 1:(length(r.common)-1)){
  chr1 <- as.character(LocMatrix[i,1])
  pos1 <- as.numeric(LocMatrix[i,2])
  
  index <- (i+1):(length(r.common))
  chr2 <- as.character(LocMatrix[index,1])
  pos2 <- as.numeric(LocMatrix[index,2])
  
  distance <- abs(pos1 - pos2)
  
  index.target <- which(chr1 == chr2 & distance > MinDis & distance < MaxDis)
  
  map1.intra <- c(map1.intra, as.numeric(map1[i, index.target]))
  map2.intra <- c(map2.intra, as.numeric(map2[i, index.target]))
  
}

t.intra <- cor(map1.intra, map2.intra, use = "pairwise.complete.obs")



map1.inter <- c()
map2.inter <- c()
for(i in 1:length(Chromosomes)){
  CHR1 <- Chromosomes[i]
  m1 <- which(LocMatrix[,1] == CHR1)
  for(j in (i+1):length(Chromosomes)){
    if(i+1 > length(Chromosomes)){
      next
    }
    CHR2 <- Chromosomes[j]
    m2 <- which(LocMatrix[,1] == CHR2)
    map1.inter <- c(map1.inter, as.numeric(map1[m1,m2]))
    map2.inter <- c(map2.inter, as.numeric(map2[m1,m2]))
  }
}
t.inter <- cor(map1.inter, map2.inter, use = "pairwise.complete.obs")

cat("Total intra data:", length(map1.intra), "\n", sep="\t")
cat("ALL:", t.all, "Intra:", t.intra, "Inter:", t.inter, "\n", sep="\t")


if(eval(parse(text=as.character(opt["out_intra"])))){
  data.x <- map1.intra
  data.y <- map2.intra
}else{
  data.x <- map1.inter
  data.y <- map2.inter
}



FILE_OUT <- as.character(opt["out"])
if(as.character(opt["height"])=="NA"){
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    h <- 5
  }else{
    h <- 500
  }
}else{
  h <- as.numeric(as.character(opt["height"]))
}
if(as.character(opt["width"])=="NA"){
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    w <- 5
  }else{
    w <- 500
  }
}else{
  w <- as.numeric(as.character(opt["width"]))
}
if(FILE_OUT != "NA"){
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    postscript(file=FILE_OUT, horizontal=FALSE, onefile=FALSE, paper="special", height=h, width=w, family="Helvetica")
  }
  if(sum((grep("\\.png$", FILE_OUT))) == 1){
    png(filename=FILE_OUT, width=w, height=h)
  }
  
  if(eval(parse(text=as.character(opt["fill"])))){
    par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  }else{
    par(oma=c(0,0,0,0), mar=c(5,5,1,1))
  }
  
  pp <- 'p'
  if(eval(parse(text=as.character(opt["NoPlot"])))){
    pp <- 'n'
  }
  plot(data.x, data.y, pch=20, type=pp, log="xy",
       xlab=as.character(opt["title1"]), ylab=as.character(opt["title2"]), cex=0.5, cex.axis=1.5, cex.lab=2.0)
  garbage <- dev.off()
}




