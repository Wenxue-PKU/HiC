# Masking TAD

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="output file from Draw_borderStrength.R"),
  make_option(c("-o", "--out"), default="NA", help="output .rds file "),
  make_option(c("--image"), default="NA", help="output mask location as image")
)
opt <- parse_args(OptionParser(option_list=option_list))

# test data
#FILE_TAD <- "W:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-04-12_domain_statistics/time2_r1/Big_ALL.txt"

options(scipen=10)
FILE_TAD <- as.character(opt["in"])
TAD <- as.matrix(read.table(FILE_TAD, header=FALSE, check.names = FALSE))
colnames(TAD) <- c("chr", "start", "end", "borderStrength", "Normalized_borderStrength", "bd", "id", "tad")
rownames(TAD) <- paste(as.character(TAD[,1]), as.numeric(TAD[,2]), as.numeric(TAD[,3]), sep=":")

IntraTAD <- function(x){
  score <- rep(FALSE, NUM_ROW)
  if(as.numeric(TAD[x, "tad"]) == 1){
    con <- ( as.numeric(TAD[, "tad"]) == 1 & (as.numeric(TAD[x,"id"]) == as.numeric(TAD[, "id"])) & (as.character(TAD[x,"chr"]) == as.character(TAD[, "chr"])) )
    score[con] <- TRUE
  }
  score
}

NUM_ROW <- nrow(TAD)
map <- cbind(sapply(1:NUM_ROW, IntraTAD))
rownames(map) <- rownames(TAD)
colnames(map) <- rownames(TAD)

FILE_object <- as.character(opt["out"])
if(FILE_object != "NA"){
  saveRDS(map, FILE_object)
}


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
  half <- upper.tri(map, diag=TRUE)
  map[half] <- NA
  
  png(file=FILE_image, width=1000, height=1000, units="px", bg="transparent")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  image(Conv02NA(Transform(map)), col=rgb(0,1,0,0.4), axes=F)
  dummy <- dev.off()
}



