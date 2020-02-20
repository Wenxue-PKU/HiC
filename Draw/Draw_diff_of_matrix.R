# diff of contact matrixを計算する
suppressWarnings(suppressMessages(library("spatstat")))

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--file1"),help="matrix file1"),
  make_option(c("-b", "--file2"),help="matrix file2"),
  make_option(c("-o", "--out"),default="NULL",help="output png file"),
  make_option(c("--matrix"), default="NULL", help="output matrix"),
  make_option(c("--adjust"), default=TRUE, help="adjust total read of file2 to file1"),
  make_option(c("--extract_first"), default=FALSE, help="In default, entire chromosome were considered for adjusting. With this option, Extract target area first then normalize"),
  make_option(c("--normalize"), default=TRUE, help="normaliz map score"),
  make_option(c("--blur"), default=TRUE, help="output png with blur"),
  make_option(c("--sigma"), default=".5", help="sigma value of Blur"),
  make_option(c("--chr"), default="all", help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--chr2"), default="NULL", help="chromosome2 name. defula is same to chr"),
  make_option(c("--start2"), default="NULL", help="start2 position. default is same to start"),
  make_option(c("--end2"), default="NULL", help="end2 position. default is same to end"),
  make_option(c("--unit"), default="v", help="unit to define score threshold p:percent or v:value"),
  make_option(c("-t", "--threshold"), default="NULL", help="min and max threshold"),
  make_option(c("--method"), default="diff", help="calculation method. diff or ratio"),
  make_option(c("--width"), default="1000", help="width of output figure"),
  make_option(c("--height"), default="NULL", help="height of output figure"),
  make_option(c("--linev_chr"), default="NULL", help="location of vertical line , separated"),
  make_option(c("--linev_pos"), default="NULL", help="location of vertical line , separated"),
  make_option(c("--lineh_chr"), default="NULL", help="location of horizontal line , separated"),
  make_option(c("--lineh_pos"), default="NULL", help="location of horizontal line , separated"),
  make_option(c("--color"), default="gentle", help="gentle or blight"),
  make_option(c("--circle"), default="NULL", help="location pairs to draw circles on output")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_matrix1 <- as.character(opt["file1"])
FILE_matrix2 <- as.character(opt["file2"])
FILE_object1 <- sub("matrix", "rds", FILE_matrix1)
FILE_object2 <- sub("matrix", "rds", FILE_matrix2)

#=========================================================
# load map1
#=========================================================
if(!file.exists(FILE_object1)){
  map1 <- as.matrix(read.table(FILE_matrix1, header=TRUE, check.names = FALSE))
}else{
  map1 <- readRDS(FILE_object1)
}
r1 <- rownames(map1)
LocList <- strsplit(r1, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

map1 <- ifelse(is.infinite(map1), NA, map1)


#=========================================================
# region check
#=========================================================
CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(CHR=="all"){
  chromosomes <- unique(as.character(LocMatrix[,1]))
  chromsome_length <- rep(0, length(chromosomes))
  names(chromsome_length) <- chromosomes
  for(c in chromosomes){
    chromsome_length[c] <- max(as.numeric(LocMatrix[as.character(LocMatrix[,1]) == c,3]))
  }
  chromosomes.sort <- sort(chromsome_length, decreasing = TRUE)
  Region <- c()
  LINE_for_chromosome_border <- c()
  for(c in names(chromosomes.sort)){
    Region <- c(Region, r1[as.character(LocMatrix[,1]) == c])
    LINE_for_chromosome_border <- c(LINE_for_chromosome_border, sum(as.character(LocMatrix[,1]) == c))
  }
  Region2 <- Region
}else{
  if(as.character(opt["end"]) == "all"){
    SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
    END <- max(as.numeric(LocMatrix[SameChromosome,3]));
  }else{
    END <- as.numeric(as.character(opt["end"]))
  }
  Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,2]) >= START) & (as.numeric(LocMatrix[,3]) <= END))
  
  if(as.character(opt["chr2"]) != "NULL"){
    CHR2 <- as.character(opt["chr2"])
  }else{
    CHR2 <- CHR
  }
  if(as.character(opt["start2"]) != "NULL"){
    START2 <- as.numeric(as.character(opt["start2"]))
  }else{
    START2 <- START
  }
  if(as.character(opt["end2"]) == "all"){
    SameChromosome <- which(as.character(LocMatrix[,1]) == CHR2)
    END2 <- max(as.numeric(LocMatrix[SameChromosome,3]));
  }else if(as.character(opt["end2"]) != "NULL"){
    END2 <- as.numeric(as.character(opt["end2"]))
  }else{
    END2 <- END
  }
  Region2 <- which((as.character(LocMatrix[,1]) == CHR2) & (as.numeric(LocMatrix[,2]) >= START2) & (as.numeric(LocMatrix[,3]) <= END2))
}


#=========================================================
# load map2
#=========================================================
if(!file.exists(FILE_object2)){
  map2 <- as.matrix(read.table(FILE_matrix2, header=TRUE, check.names = FALSE))
}else{
  map2 <- readRDS(FILE_object2)
}
map2 <- ifelse(is.infinite(map2), NA, map2)


### 先に領域を絞る
if(eval(parse(text=opt["extract_first"]))){
  map1 <- map1[Region, Region2]
  r2 <- rownames(map2)
  Region <- intersect(rownames(map1), r2)
  Region2 <- intersect(colnames(map1), r2)
  map2 <- map2[Region, Region2]
  map1 <- map1[Region, Region2]
}


# map2の合計値をmap1の合計値と同じになるように調整する
if(eval(parse(text=opt["adjust"]))){
  TotalRead <- sum(map1, na.rm=TRUE)
  map2 <- map2 / sum(map2, na.rm=TRUE) * TotalRead
}

# map1とmap2それぞれ、平均値を１にする
if(eval(parse(text=opt["normalize"]))){
  d <- nrow(map1)
  map1 <- map1 / sum(map1, na.rm=TRUE) * d * d
  map2 <- map2 / sum(map2, na.rm=TRUE) * d * d
}


# 領域を絞る（通常の順序)
if(!eval(parse(text=opt["extract_first"]))){
  map1 <- map1[Region, Region2]
  r2 <- rownames(map2)
  Region <- intersect(rownames(map1), r2)
  Region2 <- intersect(colnames(map1), r2)
  map2 <- map2[Region, Region2]
  map1 <- map1[Region, Region2]
}


#=========================================================
# 2つのmapの違いを計算する
#=========================================================
if(as.character(opt["method"]) == "diff"){
  mat.diff <- map1 - map2
}else if(as.character(opt["method"]) == "ratio"){
  mat.diff <- log2(map1 / map2)
  mat.diff <- ifelse(is.infinite(mat.diff), NA, mat.diff)
}


#=========================================================
# output matrix
#=========================================================
if(as.character(opt["matrix"]) != "NULL"){
  write.table(mat.diff, file=as.character(opt["matrix"]), quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA)
}
if(as.character(opt["out"]) == "NULL"){
  q()
}

#=========================================================
# blur image
#=========================================================
if(eval(parse(text=opt["blur"]))){
  sigma <- as.numeric(as.character(opt["sigma"]))
  t <- blur(as.im(mat.diff), sigma=sigma, bleed=FALSE)
  mat.diff <- t$v
  rownames(mat.diff) <- Region
  colnames(mat.diff) <- Region2
}


#=========================================================
# 閾値で値を区切る
#=========================================================
if(as.character(opt["unit"]) == "p"){
  NUM <- sort(abs(as.numeric(mat.diff)), na.last = NA)
  if(as.character(opt["threshold"]) == "NULL"){
    Threshold <- max(NUM, na.rm=TRUE)
  }else{
    pct <- as.numeric(as.character(opt["threshold"]))
    if(pct > 1){
      pct <- pct / 100
    }
    rank <- round(length(NUM)*pct)
    if(rank == 0){
      rank = rank +1
    }
    Threshold <- NUM[rank]
  }
}else if(as.character(opt["unit"]) == "v"){
  if(as.character(opt["threshold"]) == "NULL"){
    Threshold <- max(abs(mat.diff), na.rm=TRUE)
  }else{
    Threshold <- as.numeric(as.character(opt["threshold"]))
  }
}

mat.diff <- ifelse(mat.diff < (Threshold * -1), Threshold * -1, mat.diff)
mat.diff <- ifelse(mat.diff > Threshold, Threshold, mat.diff)


tmp <- as.numeric(mat.diff)
tmp <- tmp[!is.na(tmp)]
bk <- seq(-Threshold, Threshold, length.out=100)
mat.diff <- matrix(as.integer(cut(mat.diff, breaks = bk, include.lowest = TRUE)), nrow = nrow(mat.diff), ncol=ncol(mat.diff))

# 色分け
suppressWarnings(suppressMessages(library("RColorBrewer")))
if(as.character(opt["color"]) == "gentle"){
  pallete <- rev(brewer.pal(10, "RdBu"))  # gentle color
}
if(as.character(opt["color"]) == "blight"){
  pallete <- c(rev(brewer.pal(9, "Blues"))[-9], brewer.pal(9, "Reds"))
}
colors <- colorRampPalette(pallete)(length(bk))
colors <- colors[min(mat.diff, na.rm=TRUE):max(mat.diff, na.rm=TRUE)]


width <- as.numeric(as.character(opt["width"]))
if(as.character(opt["height"]) == "NULL"){
  height <- width / ncol(mat.diff) * nrow(mat.diff)
}else{
  height <- as.numeric(as.character(opt["height"]))
}


Transform <- function(mat){
  d = dim(mat)[1]
  n = dim(mat)[2]
  mat.rev = t(mat[d+1-c(1:d), ])
  mat.rev
}


if(as.character(opt["out"]) != "NULL"){
  png(file=as.character(opt["out"]), width=width, height=height, units="px", bg="white")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  image(Transform(mat.diff), col=colors, axes=F)
  
  if(as.character(opt["circle"]) != "NULL"){
    DATA_circle <- read.table(as.character(opt["circle"]), header=F, sep="\t", check.names = F)
    for(i in 1:nrow(DATA_circle)){
      if(DATA_circle[i,1] %in% Region & DATA_circle[i,2] %in% Region2){
        par(new=T)
        plot(which(DATA_circle[i,1] == Region), nrow(mat.diff) - which(DATA_circle[i,2] == Region2)+1, pch=21, xlim=c(1,ncol(mat.diff)), 
             ylim=c(1,nrow(mat.diff)), xaxs="i", yaxs="i", cex=2, axes=F, col='black', lwd=2.5)
      }
    }
  }
  
  if(as.character(opt["linev_chr"])!="NULL"){
    for(M in unlist(strsplit(as.character(opt["linev_pos"]), ","))){
      line <- as.numeric(gsub(" ", "", M, fixed = TRUE))
      L1 <- which((LocMatrix[,1] == as.character(opt["linev_chr"])) & (as.numeric(LocMatrix[,2]) <= line) & (as.numeric(LocMatrix[,3]) >= line))
      target <- (L1 - min(Region)) / (length(Region)-1)
      abline(h=(1-target), col="chartreuse4", lty=2, lwd=5)
    }
  }
  if(as.character(opt["lineh_chr"])!="NULL"){
    for(M in unlist(strsplit(as.character(opt["lineh_pos"]), ","))){
      line <- as.numeric(gsub(" ", "", M, fixed = TRUE))
      L2 <- which((LocMatrix[,1] == as.character(opt["lineh_chr"])) & (as.numeric(LocMatrix[,2]) <= line) & (as.numeric(LocMatrix[,3]) >= line))
      target <- (L2 - min(Region2)) / (length(Region2)-1)
      abline(v=target, col="sienna4", lty=2, lwd=5)
    }
  }
  
  if(CHR=="all"){
    Location <- cumsum(LINE_for_chromosome_border[1:(length(LINE_for_chromosome_border)-1)])/nrow(mat.diff)
    abline(v=Location, col="black", lty=1, lwd=2)
    abline(h=1-Location, col="black", lty=1, lwd=2)
  }
  
  dummy <- dev.off()
}



