# 各距離についてポアソン分布を仮定して、確率がxxx% 以下であるペアのみのマスクを出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output mask object"),
  make_option(c("-p", "--pval"), default="0.01", help="p-value to extract"),
  make_option(c("--image"), default="NA", help="output mask location as image")
)
opt <- parse_args(OptionParser(option_list=option_list))


pval <- as.numeric(as.character(opt["pval"]))

if(pval > 1){
  cat("pval should be smaller than 1\n")
  q()
}

FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}

r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
NUM_LINE <- nrow(map)


# interchrは除く
intraChr <- function(x){
  as.character(LocMatrix[x,1]) == as.character(LocMatrix[,1])
}
mask.intraChr <- cbind(sapply(1:NUM_LINE, intraChr))
rownames(mask.intraChr) <- r
colnames(mask.intraChr) <- r
map[!mask.intraChr] <- NA

# 1行目のintraChrの数を調べて、最大の距離の目安にする
NUM_maximum_intra_LINE <- sum(mask.intraChr[1,])


mask.allFalse <- mask.intraChr
mask.allFalse[,] <- FALSE

# 出力用のmatrix
mask.significant <- mask.allFalse

# 各距離で上位xx%のみTRUE、それ以外はFALSEにする
for(d in 1:(NUM_maximum_intra_LINE-1)){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  index4 <- cbind(index2, index1)
  index34 <- rbind(index3, index4)
  Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
  
  mask.temp <- mask.allFalse
  mask.temp[index34] <- TRUE

  # 指定したPvalよりも小さい値を示すためのしきい値を計算する
  Threshold <- qpois(1-pval, lambda=Average)
  
  # しきい値よりも高い部分をTRUEにする
  mask.significant[map > Threshold & mask.temp] <- TRUE
}


FILE_out <- as.character(opt["out"])
saveRDS(mask.significant, FILE_out)



Transform <- function(mat){
  d = dim(mat)[1]
  n = dim(mat)[2]
  mat.rev = t(mat[d+1-c(1:d), ])
  mat.rev
}

Conv02NA <- function(mat){
  ifelse(mat == TRUE, 1, NA)
}

FILE_image <- as.character(opt["image"])
if(FILE_image != "NA"){
  half <- upper.tri(mask.significant, diag=TRUE)
  mask.significant[half] <- NA
  
  png(file=FILE_image, width=1000, height=1000, units="px", bg="transparent")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  image(Conv02NA(Transform(mask.significant)), col=rgb(0,1,0,0.4), axes=F)
  dummy <- dev.off()
}
