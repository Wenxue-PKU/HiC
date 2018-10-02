# matrixからdistanceとスコアの関係を計算する


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output png file"),
  make_option(c("--text"), default="NA", help="output values in text(s)"),
  make_option(c("--normalize"), default="Probability", help="Average will be 1, Probability: score were divided by total read, NA:without normalization"),
  make_option(c("--mask"), default="NA", help="mask .rds file"),
  make_option(c("--ratio"), default="FALSE", help="plotting log2 ratio from first data"),
  make_option(c("--linearX"), default="FALSE", help="linear X-scale instead of log x"),
  make_option(c("--width"), default="NA", help="width of output figure"),
  make_option(c("--height"), default="NA", help="height of output figure"),
  make_option(c("--xmin"), default="NA", help="minimum distance"),
  make_option(c("--xmax"), default="NA", help="maximum distance"),
  make_option(c("--ymin"), default="NA", help="minimum score"),
  make_option(c("--ymax"), default="NA", help="maximum score"),
  make_option(c("--name"), default="NA", help="name lists"),
  make_option(c("--color"), default="NA", help="color lists if not default"),
  make_option(c("--fill"), default="FALSE", help="TRUE if not axis")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_mask <- as.character(opt["mask"])
if(FILE_mask != "NA"){
  mask <- readRDS(FILE_mask)
}

if(as.character(opt["color"]) == "NA"){
  c <- c("#e6194b", "#3cb44b", "#ffe119", "#0082c8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#fabebe", "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000080", "#808080", "#FFFFFF", "#000000")
}else{
  c <- unlist(strsplit(as.character(opt["color"]), ","))
}


getDistanceCurve <- function(FILE_matrix){
  FILE_object <- sub(".matrix", ".rds", FILE_matrix)
  if(!file.exists(FILE_object)){
    map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  }else{
    map <- readRDS(FILE_object)
  }

  # Normalization
  Normalization = as.character(opt["normalize"])
  if(Normalization == "Average"){
    d <- nrow(map)
    map <- map / sum(map, na.rm=TRUE) * d * d
  }else if(Normalization == "Probability"){
    map <- map / sum(map, na.rm=TRUE)
  }else if(Normalization == "NA"){
    # don't change anything
  }else{
    cat("Normalization parameter is wrong")
    q()
  }
  
  if(FILE_mask != "NA"){
    r.common <- intersect(rownames(map), rownames(mask))
    map <- map[r.common, r.common]
    mask <- mask[r.common, r.common]
    map[!mask] <- NA
  }
  
  LocList <- strsplit(rownames(map), ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
  NUM_LINE <- nrow(map)
  
  # inter-chromosomeは除外する
  IntraChr <- function(x){
    as.character(LocMatrix[x,1]) == as.character(LocMatrix[,1])
  }
  mask_intra <- cbind(sapply(1:NUM_LINE, IntraChr))
  map[!mask_intra] <- NA
  

  x <- c()
  y <- c()
  
  for(d in 0:(NUM_LINE-1)){
    index1 <- 1:(NUM_LINE - d)
    index2 <- index1 + d
    index3 <- cbind(index1, index2)
    # index4 <- cbind(index2, index1)
    Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
    if(!is.na(Average) & Average != 0){
      distance <- d*Resolution
      if(distance == 0){
        distance = 1
      }
      x <- c(x, distance)
      y <- c(y, Average)
    }
    # cat(Average, "\n")
  }
  Data <- list()
  Data[["x"]] <- x
  Data[["y"]] <- y
  Data
}


FILEs <- unlist(strsplit(as.character(opt["in"]), ","))
if(as.character(opt["text"]) != "NA"){
  FILE_text <- unlist(strsplit(as.character(opt["text"]), ","))
}
SAMPLE_NUMBER <- length(FILEs)
Dis_DATA <- list()
for(i in 1:SAMPLE_NUMBER){
  Dis_DATA[[i]] <- getDistanceCurve(FILEs[i])
  
  # output
  if(as.character(opt["text"]) != "NA"){
    OUT_TEXT <- cbind(Dis_DATA[[i]][["x"]], Dis_DATA[[i]][["y"]])
    write.table(OUT_TEXT, FILE_text[i], row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE, append = FALSE, eol = "\n")
  }
}

if(as.character(opt["out"]) == "NA"){
  q()
}


# 1つめのファイルとのlog2 ratioを計算する
if(eval(parse(text=as.character(opt["ratio"])))){
  Dis_DATA_new <- list()
  X_in_1st_graph <- as.numeric(Dis_DATA[[1]][["x"]])
  Y_in_1st_graph <- as.numeric(Dis_DATA[[1]][["y"]])
  
  for(i in 2:SAMPLE_NUMBER){
    X_in_ith_graph <- as.numeric(Dis_DATA[[i]][["x"]])
    Y_in_ith_graph <- as.numeric(Dis_DATA[[i]][["y"]])
    x <- c()
    y <- c()
    for(k in 1:length(X_in_ith_graph)){
      if(X_in_ith_graph[k] %in% X_in_1st_graph){
        index <- which(X_in_1st_graph == X_in_ith_graph[k])
        x <- c(x, X_in_ith_graph[k])
        y <- c(y, log2(Y_in_ith_graph[k] / Y_in_1st_graph[index]))
      }
    }
    Data <- list()
    Data[["x"]] <- x
    Data[["y"]] <- y
    Dis_DATA_new[[i-1]] <- Data
  }
  Dis_DATA <- Dis_DATA_new
  SAMPLE_NUMBER <- length(FILEs) -1
  rm(Dis_DATA_new, Data)
}


FILE_OUT <- as.character(opt["out"])
if(as.character(opt["height"])=="NA"){
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    h <- 4
  }else{
    h <- 12
  }
}else{
  h <- as.numeric(as.character(opt["height"]))
}
if(as.character(opt["width"])=="NA"){
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    w <- 4
  }else{
    w <- 12
  }
}else{
  w <- as.numeric(as.character(opt["width"]))
}





if(as.character(opt["xmin"]) == "NA"){
  xmin <- min(as.numeric(Dis_DATA[[1]][["x"]]), na.rm = TRUE)
}else{
  xmin <- as.numeric(as.character(opt["xmin"]))
}

if(as.character(opt["xmax"]) == "NA"){
  xmax <- max(as.numeric(Dis_DATA[[1]][["x"]]), na.rm = TRUE)
}else{
  xmax <- as.numeric(as.character(opt["xmax"]))
}

if(as.character(opt["ymin"]) == "NA"){
  ymin <- min(as.numeric(Dis_DATA[[1]][["y"]]), na.rm = TRUE) * 0.8
}else{
  ymin <- as.numeric(as.character(opt["ymin"]))
}

if(as.character(opt["ymax"]) == "NA"){
  ymax <- max(as.numeric(Dis_DATA[[1]][["y"]]), na.rm = TRUE) * 1.2
}else{
  ymax <- as.numeric(as.character(opt["ymax"]))
}



if(sum((grep("\\.eps$", FILE_OUT))) == 1){
  postscript(file=FILE_OUT, horizontal=FALSE, onefile=FALSE, paper="special", height=h, width=w, family="Helvetica")
  cex.axis=0.7
  cex.lab=0.8
}
if(sum((grep("\\.png$", FILE_OUT))) == 1){
  png(filename=FILE_OUT, width=w, height=h, bg="white", res = 100, units = "cm")
  cex.axis=1.4
  cex.lab=1.5
}


LogAxis <- "xy"
Ylabels <- as.character(opt["normalize"])
if(eval(parse(text=as.character(opt["linearX"])))){
  LogAxis <- "y"
}
if(eval(parse(text=as.character(opt["ratio"])))){
  LogAxis <- "x"
  Ylabels <- "Log2 ratio of probability"
  if(eval(parse(text=as.character(opt["linearX"])))){
    LogAxis <- ""
  }
}


if(eval(parse(text=as.character(opt["fill"])))){
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  plot(Dis_DATA[[1]][["x"]], Dis_DATA[[1]][["y"]], type='l', xaxs = "i", yaxs="i", xlab="Distance", ylab=Ylabels, 
       col=c[1], log = LogAxis, cex.axis=cex.axis, cex.lab=cex.lab, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), bty = "n")
  if(SAMPLE_NUMBER > 1){
    for(n in 2:SAMPLE_NUMBER){
      par(new=T)
      plot(Dis_DATA[[n]][["x"]], Dis_DATA[[n]][["y"]], type='l', xaxs = "i", yaxs="i", axes=F, xlab="", ylab="", col=c[n], 
           log = LogAxis, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax), bty = "n")
    }
  }
}else{
  par(oma=c(0,0,0,0), mar=c(5,5,1,1))
  plot(Dis_DATA[[1]][["x"]], Dis_DATA[[1]][["y"]], type='l', xaxs = "i", yaxs="i", xlab="Distance", ylab=Ylabels,
       col=c[1], log = LogAxis, cex.axis=cex.axis, cex.lab=cex.lab, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax))
  if(SAMPLE_NUMBER > 1){
    for(n in 2:SAMPLE_NUMBER){
      par(new=T)
      plot(Dis_DATA[[n]][["x"]], Dis_DATA[[n]][["y"]], type='l', xaxs = "i", yaxs="i", axes=F, xlab="", ylab="", col=c[n], 
           log = LogAxis, lwd=2, xlim=c(xmin, xmax), ylim=c(ymin, ymax))
    }
  }
}
if(as.character(opt["name"]) != "NA"){
  names <- unlist(strsplit(as.character(opt["name"]), ","))
  legend("bottomleft", legend=names, lwd=2, lty=1, col=c[1:length(names)], bty='n')
}
if(eval(parse(text=as.character(opt["ratio"])))){
  abline(h=0, col='blue', lty=3, lwd=2)
}
grabage <- dev.off()






