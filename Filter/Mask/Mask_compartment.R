# Masking Compartment

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="output location file from PCA_analysis.R"),
  make_option(c("-o", "--out"), default="NA", help="output .rds file "),
  make_option(c("--image"), default="NA", help="output mask location as image")
)
opt <- parse_args(OptionParser(option_list=option_list))


# test data
#FILE_DATA <- "W:/Project/2015-10-13_Human_youngOld_IMR90/out/2017-05-10_summary_of_current_domain/PCA/location/SCIMR_chr1.txt"

options(scipen=10)
FILE_DATA <- as.character(opt["in"])
DATA <- as.matrix(read.table(FILE_DATA, header=FALSE, check.names = FALSE))
colnames(DATA) <- c("chr", "start", "end", "pca_score", "Compartment", "bd")
rownames(DATA) <- paste(as.character(DATA[,1]), as.numeric(DATA[,2]), as.numeric(DATA[,3]), sep=":")

IntraCompartment <- function(x){
  score <- rep(NA, NUM_ROW)
  targetCompartment <- as.character(DATA[x,"Compartment"])
  con <- ( targetCompartment == as.character(DATA[,"Compartment"]) & (as.character(DATA[x,"chr"]) == as.character(DATA[, "chr"])) )
  score[con] <- targetCompartment
  score
}


NUM_ROW <- nrow(DATA)
map <- cbind(sapply(1:NUM_ROW, IntraCompartment))
rownames(map) <- rownames(DATA)
colnames(map) <- rownames(DATA)

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

Conv2score <- function(mat){
  mat <- ifelse(mat == "A", NA, mat)
  mat <- ifelse(mat == "B", 1, mat)
  mat
}

FILE_image <- as.character(opt["image"])
if(FILE_image != "NA"){
  half <- upper.tri(map, diag=TRUE)
  map[half] <- NA
  map.draw <- Conv2score(Transform(map))
  mode(map.draw) <- 'numeric'
  
  png(file=FILE_image, width=1000, height=1000, units="px", bg="transparent")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  image(map.draw, col=rgb(0,1,0,0.4), axes=F)
  dummy <- dev.off()
}