# Draw mask region

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help=".rds file of mask"),
  make_option(c("--triangle"), default="TRUE", help="plot only half of triangle"),
  make_option(c("--chr"),help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--color"), default="NA", help="color"),
  make_option(c("-o", "--out"), help="output mask location as image")
)
opt <- parse_args(OptionParser(option_list=option_list))

options(scipen=10)
FILE_mask <- as.character(opt["in"])
mask <- readRDS(FILE_mask)
r <- rownames(mask)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

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

CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))

if(CHR=="all"){
  mask <- SeparatePombe(mask, LocMatrix)
}else{
  if(as.character(opt["end"]) == "all"){
    SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
    END <- max(as.numeric(LocMatrix[SameChromosome,3]));
  }else{
    END <- as.numeric(as.character(opt["end"]))
  }
  Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,3]) >= START) & (as.numeric(LocMatrix[,2]) <= END))
  mask <- mask[Region, Region]
}




if(eval(parse(text=as.character(opt["triangle"])))){
  half <- upper.tri(mask, diag=TRUE)
  mask[half] <- NA
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

COLOR <- as.character(opt["color"])
if(COLOR == "NA"){
  COLOR=rgb(0,1,0,0.4)
}

FILE_image <- as.character(opt["out"])
png(file=FILE_image, width=1000, height=1000, units="px", bg="transparent")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
image(Conv02NA(Transform(mask)), col=COLOR, axes=F)
dummy <- dev.off()




