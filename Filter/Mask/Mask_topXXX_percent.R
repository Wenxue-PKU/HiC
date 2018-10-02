# 上位xxx%を出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output mask object"),
  make_option(c("-p", "--percent"), default="0.01", help="parcent to extract value should be 0 - 1"),
  make_option(c("--image"), default="NA", help="output mask location as image")
)
opt <- parse_args(OptionParser(option_list=option_list))


percent <- as.numeric(as.character(opt["percent"]))

if(percent > 1){
  cat("percent should be smaller than 1\n")
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

# 出力用のmatrix
mask.significant<- map
mode(mask.significant) <- "logical"
mask.significant[,] <- FALSE

scores <- as.numeric(map)
scores.sort <- sort(scores, decreasing = TRUE, na.last = NA)
Threshold <- scores.sort[as.integer(length(scores.sort) * percent)]

mask.significant[map > Threshold] <- TRUE


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
