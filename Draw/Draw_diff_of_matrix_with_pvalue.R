# diff of contact matrixを計算する
#######################################
# file1 の方が強い場合赤色になる
#######################################



# .libPaths(.libPaths()[c(2,3,1)])
suppressWarnings(suppressMessages(library("spatstat")))

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--file1"),help="matrix file1"),
  make_option(c("-b", "--file2"),help="matrix file2"),
  make_option(c("--bio1"),help="matrix file biological replica1"),
  make_option(c("--bio2"),help="matrix file biological replica2"),
  make_option(c("--adjust_control"), default=FALSE, help="half of the ratio were multiplied by -1"),
  make_option(c("--FDR"), default="0.05", help="FDR threshold"),
  make_option(c("-o", "--out"),default="NULL",help="output png file"),
  make_option(c("--matrix"), default="NULL", help="output matrix"),
  make_option(c("--chr"), default="all", help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--chr2"), default="NULL", help="chromosome2 name. defula is same to chr"),
  make_option(c("--start2"), default="NULL", help="start2 position. default is same to start"),
  make_option(c("--end2"), default="NULL", help="end2 position. default is same to end"),
  make_option(c("--width"), default="1000", help="width of output figure"),
  make_option(c("--height"), default="NULL", help="height of output figure"),
  make_option(c("--hist"), default="NA", help="output file name of pvalue histogram"),
  make_option(c("--linev_chr"), default="NULL", help="location of vertical line , separated"),
  make_option(c("--linev_pos"), default="NULL", help="location of vertical line , separated"),
  make_option(c("--lineh_chr"), default="NULL", help="location of horizontal line , separated"),
  make_option(c("--lineh_pos"), default="NULL", help="location of horizontal line , separated")
)
opt <- parse_args(OptionParser(option_list=option_list))

TakeMiddleV <- function(mat, minVal, maxVal){
  mat.new <- ifelse(mat < minVal, minVal, mat)
  ifelse(mat.new > maxVal, maxVal, mat.new)
}

Transform <- function(mat){
  d = dim(mat)[1]
  n = dim(mat)[2]
  mat.rev = t(mat[d+1-c(1:d), ])
  mat.rev
}

SeparatePombe <- function(m, L){
  blankNumber = as.integer(dim(m)[1]*0.01)  
  blank <- rep(NA, dim(m)[1])
  
  chr1 <- which(as.character(L[,1]) == "I")
  chr2 <- which(as.character(L[,1]) == "II")
  chr3 <- which(as.character(L[,1]) == "III")
  
  m.rbind <- rbind(m[chr1,], matrix(rep(blank,blankNumber), nrow=blankNumber), m[chr2,], matrix(rep(blank,blankNumber),nrow=blankNumber), m[chr3,])
  m.rbind
  rowNumber <- dim(m.rbind)[1]
  blank <- rep(NA, rowNumber)
  m.cbind <- cbind(m.rbind[,chr1], matrix(rep(blank,blankNumber), ncol=blankNumber), m.rbind[,chr2], matrix(rep(blank,blankNumber), ncol=blankNumber), m.rbind[,chr3])
  m.cbind
}

SeparatePombe_name <- function(m, L){
  blank <- rep(NA, length(m) * 0.01)
  
  chr1 <- which(as.character(L[,1]) == "I")
  chr2 <- which(as.character(L[,1]) == "II")
  chr3 <- which(as.character(L[,1]) == "III")
  
  c(m[chr1], blank, m[chr2], blank, m[chr3])
}


FILE_matrix1 <- as.character(opt["file1"])
FILE_matrix2 <- as.character(opt["file2"])
FILE_matrix3 <- as.character(opt["bio1"])
FILE_matrix4 <- as.character(opt["bio2"])
FILE_object1 <- sub("matrix", "rds", FILE_matrix1)
FILE_object2 <- sub("matrix", "rds", FILE_matrix2)
FILE_object3 <- sub("matrix", "rds", FILE_matrix3)
FILE_object4 <- sub("matrix", "rds", FILE_matrix4)

map1 <- readRDS(FILE_object1)
map2 <- readRDS(FILE_object2)
map_bio1 <- readRDS(FILE_object3)
map_bio2 <- readRDS(FILE_object4)

r1 <- rownames(map1)
r2 <- rownames(map2)
r3 <- rownames(map_bio1)
r4 <- rownames(map_bio2)
r.common <- intersect(r1, r2)
r.common.original <- r.common
r.common2 <- intersect(r3, r4)

map1 <- map1[r.common, r.common]
map2 <- map2[r.common, r.common]
map_bio1 <- map_bio1[r.common2, r.common2]
map_bio2 <- map_bio2[r.common2, r.common2]

draw_map <- map1
draw_map[,] <- 1

LocList <- strsplit(r.common, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)


d <- length(r.common)
map1 <- map1 / sum(map1, na.rm=TRUE) * d * d
map2 <- map2 / sum(map2, na.rm=TRUE) * d * d
d <- length(r.common2)
map_bio1 <- map_bio1 / sum(map_bio1, na.rm=TRUE) * d * d
map_bio2 <- map_bio2 / sum(map_bio2, na.rm=TRUE) * d * d


CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))

target.ratio <- log2(map1/map2)
target.ratio <- ifelse(is.infinite(target.ratio), NA, target.ratio)
target.ave <- (map1 + map2)/2

control.ratio <- as.numeric(log2(map_bio1/map_bio2))
control.ave <- as.numeric((map_bio1 + map_bio2)/2)
control.strange <- is.na(control.ratio) | is.infinite(control.ratio)
control.ratio <- control.ratio[!control.strange]
control.ave <- control.ave[!control.strange]


if(eval(parse(text=as.character(opt["adjust_control"])))){
  index <- sample(length(control.ratio), length(control.ratio)/2, replace = FALSE)
  control.ratio[index] <- control.ratio[index] * -1
}



lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

# controlが最低でも1000あるようにカテゴリー分け
control.ave.sort <- sort(control.ave, na.last = NA)
index_max <- length(control.ave.sort)-5000
b <- c(0, lseq(control.ave.sort[5000], control.ave.sort[index_max], 20),1e9)

target.category <- sapply(1:nrow(target.ave), function(i){cut(target.ave[,i], breaks=b, right=FALSE)})
control.category <- cut(control.ave, breaks = b, right=FALSE)

NUM_sampling <- 10000
background <- list()
category_list <- levels(factor(control.category))
for(ca in category_list){
  background[[ca]] <- sample(control.ratio[which(control.category == ca)], NUM_sampling, replace = TRUE)
}

if(CHR=="all"){
  target.ratio <- SeparatePombe(target.ratio, LocMatrix)
  target.category <- SeparatePombe(target.category, LocMatrix)
  target.ave <- SeparatePombe(target.ave, LocMatrix)
  draw_map <- SeparatePombe(draw_map, LocMatrix)
  r.common <- SeparatePombe_name(r.common, LocMatrix)
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
  
  target.ratio <- target.ratio[Region, Region2]
  target.ave <- target.ave[Region, Region2]
  target.category <- target.category[Region, Region2]
  draw_map <- draw_map[Region, Region2]
  r.common <- r.common[Region]
}


pvalCalc <- function(cate, val){
  if(is.na(val)){
    NA
  }else if(exists(as.character(cate), where=background)){
    sum(val > background[[cate]]) / NUM_sampling
  }else{
    NA
  }
}
pval <- sapply(1:nrow(target.ratio), function(i){mapply(pvalCalc, target.category[i,], target.ratio[i,])})
pval <- ifelse(pval > 0.5, 1-pval, pval)
qval <- p.adjust(pval, method = "BH")

FILE_hist <- as.character(opt["hist"])
if(FILE_hist != "NA"){
  png(file=FILE_hist, width=400, height=300, units="px", bg="white")
  par(oma=c(4,4,1,1), mar=c(0,0,0,0))
  hist(as.numeric(pval), breaks=50, col="grey", main="")
  dummy <- dev.off()
}

FDR_threshold <- as.numeric(as.character(opt["FDR"]))
pval_threshold <- max(as.numeric(pval[qval <= FDR_threshold]), na.rm = TRUE)

draw_map<- ifelse(pval <= pval_threshold & target.ratio > 0, 3, draw_map)
draw_map <- ifelse(pval <= pval_threshold & target.ratio < 0, 2, draw_map)

# draw_map<- ifelse(pval < pval_threshold & target.ratio > 0 & target.ave > 1, 3, draw_map)
# draw_map<- ifelse(pval < pval_threshold & target.ratio > 0 & target.ave < 1, 3, draw_map)


width <- as.numeric(as.character(opt["width"]))
if(as.character(opt["height"]) == "NULL"){
  height <- width / nrow(draw_map) * ncol(draw_map)
}else{
  height <- as.numeric(as.character(opt["height"]))
}

FILE_OUT <- as.character(opt["out"])
if(FILE_OUT != "NULL"){
  png(file=FILE_OUT, width=width, height=height, units="px", bg="white")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  image(Transform(draw_map), col=c('#ffffbf', '#2c7bb6','#d7191c'), axes=F)
  
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
  dummy <- dev.off()
}


if(as.character(opt["matrix"]) != "NULL"){
  colnames(pval) <- r.common
  rownames(pval) <- r.common
  pval <- pval[r.common.original, r.common.original]
  suppressWarnings( write.table(pval, file=as.character(opt["matrix"]), quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA))
  cat("Pvale threshold to get FDR < ", FDR_threshold, " : ", pval_threshold, "\n", sep="")
}



