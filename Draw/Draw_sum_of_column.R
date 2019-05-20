# contact mapのcolumnごとのbiasを出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output png file"),
  make_option(c("--chr"),help="chromosome name all for all chromosome"),
  make_option(c("--refresh"), default="FALSE", help="FALSE: use existing map object if eixts"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--min"), default="NULL", help="minimum score for drawing"),
  make_option(c("--max"), default="NULL", help="maximu score for drawing"),
  make_option(c("--width"), default="NA", help="width of output figure"),
  make_option(c("--height"), default="NA", help="height of output figure"),
  make_option(c("--color"), default="random", help="graph color"),
  make_option(c("--fill"), default="TRUE", help="TRUE if not axis")
)
opt <- parse_args(OptionParser(option_list=option_list))

Flag_refresh <- eval(parse(text=as.character(opt["refresh"])))

FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object) || Flag_refresh){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  saveRDS(map, FILE_object)
}else{
  map <- readRDS(FILE_object)
}
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(as.character(opt["end"]) == "all"){
  END <- max(DATA[,2])
}else{
  END <- as.numeric(as.character(opt["end"]))
}

t.filter <- which(as.character(LocMatrix[,1]) == CHR & as.numeric(LocMatrix[,3]) >= START &  as.numeric(LocMatrix[,2]) <= END)

# 値をnormalizeする（各列の合計が1になるようにする）
d <- length(r)
map <- map / sum(map, na.rm=TRUE) * d

# 各columnごとの合計値を計算
sum <- apply(map, 2, sum, na.rm=TRUE)


DATA <- cbind(LocMatrix, sum)
rownames(DATA) <- r

FILE_size <-  file.info(FILE_matrix)$size
if(as.character(opt["color"]) == "random"){
  set.seed(seed = FILE_size)
  random <- as.integer(runif(1, min=1, max=100))
  c <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                          "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))( 100 )
  graphColor <- c[random]
}else{
  graphColor <- as.character(opt["color"])
}





FILE_OUT <- as.character(opt["out"])
if(as.character(opt["height"])=="NA"){
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    h <- 2.2
  }else{
    h <- 100
  }
}else{
  h <- as.numeric(as.character(opt["height"]))
}
if(as.character(opt["width"])=="NA"){
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    w <- 8
  }else{
    w <- 1000
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

DrawEach <- function(i){
  polygon(c(as.numeric(DATA[i,2]), as.numeric(DATA[i,2]), as.numeric(DATA[i,3]), as.numeric(DATA[i,3])), 
          c(as.numeric(DATA[i,4]), 0, 0, as.numeric(DATA[i,4])), col=graphColor, border=graphColor)
}

if(eval(parse(text=as.character(opt["fill"])))){
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  plot(c(START, END), c(0,1), type='n', xaxs = "i", yaxs="i", xlab="", ylab="", ylim=c(0,2), axes=FALSE, xlim=c(START, END))
}else{
  par(oma=c(0,0,0,0), mar=c(3,4,1,1))
  plot(c(START, END), c(0,1), type='n', xaxs = "i", yaxs="i", xlab="", ylab="", ylim=c(0,2), xlim=c(START, END))
}

dummy <- sapply(t.filter, DrawEach)
grabage <- dev.off()
