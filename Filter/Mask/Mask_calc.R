# Mask calculation

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--mask1"), help="mask file1"),
  make_option(c("-b", "--mask2"), help="mask file2"),
  make_option(c("-o", "--out"), help="output .rds file "),
  make_option(c("--method"), default="intersect", help="how to calculate. intersect, merge, subtract"),
  make_option(c("--image"), default="NA", help="output mask location as image")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_mask1 <- as.character(opt["mask1"])
FILE_mask2 <- as.character(opt["mask2"])
mask1 <- readRDS(FILE_mask1)
mask2 <- readRDS(FILE_mask2)
r.common <- intersect(rownames(mask1), rownames(mask2))

mask1 <- mask1[r.common, r.common]
mask2 <- mask2[r.common, r.common]


method <- as.character(opt["method"])
switch(method,
       "intersect" = map <- mask1 & mask2,
       "merge" = map <- mask1 | mask2,
       "subtract" = map <- mask1 & (!mask2),
       cat(method, " is not defined")
)


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
  png(file=FILE_image, width=nrow(map), height=nrow(map), units="px", bg="transparent")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  image(Conv02NA(Transform(map)), col=rgb(0,1,0,0.4), axes=F)
  dummy <- dev.off()
}

