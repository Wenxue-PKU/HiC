# Masking distance

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="template matrix"),
  make_option(c("-o", "--out"), help="output .rds file "),
  make_option(c("--min"), default="0", help="minimum distance"),
  make_option(c("--max"), default="999999999", help="maximum distance"),
  make_option(c("--image"), default="NA", help="output mask location as image")
)
opt <- parse_args(OptionParser(option_list=option_list))


Threshold_min <- as.numeric(as.character(opt["min"]))
Threshold_max <- as.numeric(as.character(opt["max"]))


FILE_template <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_template)
map_ref <- readRDS(FILE_object)

# falaseで初期化
map <- matrix(rep(FALSE, nrow(map_ref)*ncol(map_ref)), nrow=nrow(map_ref))
rownames(map) <- rownames(map_ref)
colnames(map) <- colnames(map_ref)


r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
NUM_LINE <- nrow(map)


for(d in 1:(NUM_LINE-1)){
  distance <- d*Resolution
  if(distance >= Threshold_min & distance <= Threshold_max){
    index1 <- 1:(NUM_LINE - d)
    index2 <- index1 + d
    index3 <- cbind(c(index1, index2), c(index2, index1))
    map[index3] <- TRUE
  }
}


# inter-chromosomeは除外する
interChr <- function(x){
  as.character(LocMatrix[x,1]) != as.character(LocMatrix[,1])
}
mask_inter <- cbind(sapply(1:NUM_LINE, interChr))
map[mask_inter] <- FALSE



FILE_object <- as.character(opt["out"])
saveRDS(map, FILE_object)


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
  png(file=FILE_image, width=NUM_LINE, height=NUM_LINE, units="px", bg="transparent")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  image(Conv02NA(Transform(map)), col=rgb(0,1,0,0.4), axes=F)
  dummy <- dev.off()
}

