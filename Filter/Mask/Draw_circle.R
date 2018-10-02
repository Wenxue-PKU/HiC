# matrixでmaskの場所を◯で囲む

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help=".rds file of mask"),
  make_option(c("-o", "--out"),help="output png file"),
  make_option(c("--chr"),help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--color"), default="white", help="color of circle"),
  make_option(c("--circle"), default="1.3", help="circle size"),
  make_option(c("--width"), default="1000", help="width of output figure"),
  make_option(c("--height"), default="NULL", help="height of output figure")
)
opt <- parse_args(OptionParser(option_list=option_list))


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

options(scipen=10)
FILE_mask <- as.character(opt["in"])
mask <- readRDS(FILE_mask)
r <- rownames(mask)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)


CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))

if(CHR=="all"){
  Region <- 1:nrow(mask)
  Region2 <- Region
  mask <- SeparatePombe(mask, LocMatrix)
}else{
  if(as.character(opt["end"]) == "all"){
    SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
    END <- max(as.numeric(LocMatrix[SameChromosome,3]));
  }else{
    END <- as.numeric(as.character(opt["end"]))
  }
  Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,3]) >= START) & (as.numeric(LocMatrix[,2]) <= END))
  
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
  Region2 <- which((as.character(LocMatrix[,1]) == CHR2) & (as.numeric(LocMatrix[,3]) >= START2) & (as.numeric(LocMatrix[,2]) <= END2))
  
  mask <- mask[Region, Region2]
}
mask <- Transform(mask)



locations <- which(mask==TRUE, arr.ind = TRUE)


unit <- 1/(dim(mask)[1])
convLoc <- function(x){
  x * unit - 0.5*unit
}

width <- as.numeric(as.character(opt["width"]))
if(as.character(opt["height"]) == "NULL"){
  height <- width / nrow(mask) * ncol(mask)
}else{
  width <- as.numeric(as.character(opt["height"]))
}


FILE_OUT <- as.character(opt["out"])
COLOR <- as.character(opt["color"])
CIRCLE_SIZE <- as.numeric(as.character(opt["circle"]))
png(file=FILE_OUT, width=width, height=height, units="px", bg="transparent")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
plot(convLoc(locations[,1]), convLoc(locations[,2]), type='p', pch=1, 
     xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), cex=CIRCLE_SIZE, axes=F, xaxs="i", yaxs="i", lwd=3, col=adjustcolor(COLOR, alpha.f = 0.7))
dummy <- dev.off()






