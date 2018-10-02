# compartment A/B を主成分分析で解析する
# plus ('blue')がopenでminus('red')がclose

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"), default="NA", help="output png file"),
  make_option(c("--quiet"), default="TRUE", help="output log or not"),
  make_option(c("--geneDensity"), default="NA", help="gene density file"),
  make_option(c("--refresh"), default="FALSE", help="FALSE: use existing map object if eixts"),
  make_option(c("--location"), default="NA", help="PCA score"),
  make_option(c("--chr"),help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--ymax"), default="NA", help="Maximum score of Y-axis"),
  make_option(c("--ymin"), default="NA", help="Minimum score of Y-axis"),
  make_option(c("--width"), default="NA", help="width of output figure"),
  make_option(c("--fill"), default="FALSE", help="fill image or not"),
  make_option(c("--height"), default="NA", help="height of output figure")
)
opt <- parse_args(OptionParser(option_list=option_list))

Flag_refresh <- eval(parse(text=as.character(opt["refresh"])))
Flag_quiet <- eval(parse(text=as.character(opt["quiet"])))

FILE_geneDensity <- as.character(opt["geneDensity"])
if(FILE_geneDensity == "NA"){
  cat("gene density infomration is required")
  q()
}


FILE_matrix <- as.character(opt["in"])
FILE_object <- sub("matrix", "rds", FILE_matrix)
if(!file.exists(FILE_object) || Flag_refresh){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  saveRDS(map, FILE_object)
}else{
  map <- readRDS(FILE_object)
}
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

RESOLUTION <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1


CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(as.character(opt["end"]) == "all"){
  END <- max(as.numeric(LocMatrix[,3]))
}else{
  END <- as.numeric(as.character(opt["end"]))
}

# Observed / Expectのmatrixに変換する
map_expect <- map
NUM_LINE <- nrow(map)
for(d in 0:(NUM_LINE-1)){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  index4 <- cbind(index2, index1)
  Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
  map_expect[index3] <- Average
  map_expect[index4] <- Average
}
map <- ifelse(map_expect == 0, 1, map / map_expect)

# もしsdが0だったらNAにする
sdlist <- apply(map,1,sd, na.rm=TRUE)
index <- which(sdlist == 0)
map[index, ] <- NA
map[, index] <- NA

getCor <- function(i){
  apply(map, 2, function(x) { cor(x, map[,i], use="pairwise.complete.obs", method="pearson")})
}
options(warn=-1)
map <- rbind(sapply(1:ncol(map), getCor))
colnames(map) <- rownames(map)
options(warn=1)

SUM <- apply(is.na(map), 2, sum)
NonNABin <- r[SUM != nrow(map)]

# PCA analysis
pca <- prcomp(map[NonNABin, NonNABin],  scale=TRUE)
if(!Flag_quiet){
  summary(pca)
}

# PCA score
pca_score <- rep(NA, nrow(map))

# 合計で何個の連続したcompartmentがあるのかを数え、もし、その数が10未満で、第2主成分では、10以上であった場合、第2主成分を用いる
FirstComponent <- pca$x[,1]
SecondComponent <- pca$x[,2]
names(pca_score) <- r

checkCompartmentNum <- function(scores){
  scores <- sign(scores)
  scores <- ifelse(is.na(scores), 0, scores)
  NUM_total_compartment <- 1
  old_sign <- scores[1]
  for(i in 2:length(scores)){
    new_sign <- scores[i]
    if(new_sign != old_sign){
      NUM_total_compartment  <- NUM_total_compartment + 1
    }
    old_sign <- new_sign
  }
  NUM_total_compartment
}

if(checkCompartmentNum(FirstComponent) < 10 && checkCompartmentNum(SecondComponent) > 10){
  pca_score[NonNABin] <- SecondComponent
}else{
  pca_score[NonNABin] <- FirstComponent
}


# 遺伝子の密度が多い方をCompartmentA、低い方をcompartmentBとする
GeneDen <- as.matrix(read.table(FILE_geneDensity, header=FALSE, check.names=FALSE, row.names=1))
group1 <- intersect(names(which(pca_score > 0)), r)
group2 <- intersect(names(which(pca_score < 0)), r)


if(mean(as.numeric(GeneDen[group1,1])) < mean(as.numeric(GeneDen[group2,1]))){
  # group2がcompartmentA
  if(sum(pca_score[group2], na.rm=TRUE) < 0){
    pca_score <- pca_score * (-1)
  }
}else{
  # group1がcompartmentA
  if(sum(pca_score[group1], na.rm=TRUE) < 0){
    pca_score <- pca_score * (-1)
  }
}


FILE_location <- as.character(opt["location"])
if(FILE_location!= "NA"){
  # AとBを定義
  Compartment <- rep("NA", length(pca_score))
  Compartment <- ifelse(pca_score > 0, "A", Compartment)
  Compartment <- ifelse(pca_score < 0, "B", Compartment)
  options(scipen=10)
  OUTPUT <- cbind(LocMatrix, pca_score,Compartment)
  write.table(OUTPUT, file=FILE_location, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE, append=FALSE)
}

FILE_OUT <- as.character(opt["out"])
if(FILE_OUT != "NA"){
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
      w <- 4
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
    png(file=FILE_OUT, width=w, height=h, units="px", bg="white")
  }
  
  NUM <- sort(pca_score, na.last = NA)
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
  
  index <- which(as.character(LocMatrix[,1]) == CHR & as.numeric(LocMatrix[,3]) >= START &  as.numeric(LocMatrix[,2]) <= END)
  x1 <- as.numeric(LocMatrix[index,2])
  x2 <- as.numeric(LocMatrix[index,3])
  y <- pca_score[index]
  
  drawPoly <- function(i){
    polygon(c(x1[i],x2[i],x2[i],x1[i]), c(y[i], y[i], 0, 0), col=ifelse(y[i]>0,'red', 'blue'), ylim=c(YMIN,YMAX), xlim=c(START,END))
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
  
  gabage <- dev.off()
}


