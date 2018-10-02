# Output correlation of two matrix

suppressPackageStartupMessages(library("hicrep"))
suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--file1"),help="matrix file1"),
  make_option(c("-b", "--file2"),help="matrix file2"),
  make_option(c("--adjust"), default=FALSE, help="adjusting read (TRUE) or not (FALSE)"),
  make_option(c("--chromosome"), default="NA", help="chromosome to caclulate"),
  make_option(c("--max_distance"), default="100000000", help="maximum distance for HiC rep anlaysis")
)
opt <- parse_args(OptionParser(option_list=option_list))




FILE_matrix1 <- as.character(opt["file1"])
FILE_matrix2 <- as.character(opt["file2"])
FILE_object1 <- sub("matrix", "rds", FILE_matrix1)
FILE_object2 <- sub("matrix", "rds", FILE_matrix2)

if(!file.exists(FILE_object1)){
  map1 <- as.matrix(read.table(FILE_matrix1, header=TRUE, check.names = FALSE))
}else{
  map1 <- readRDS(FILE_object1)
}
if(!file.exists(FILE_object2)){
  map2 <- as.matrix(read.table(FILE_matrix2, header=TRUE, check.names = FALSE))
}else{
  map2 <- readRDS(FILE_object2)
}

r <- intersect(rownames(map1), rownames(map2))
map1 <- map1[r,r]
map2 <- map2[r,r]

if(eval(parse(text=as.character(opt["adjust"])))){
  NUM1 <- sum(map1, na.rm = TRUE)
  NUM2 <- sum(map2, na.rm = TRUE)
  if(NUM1 > NUM2){
    map1 <- map1 / NUM1 * NUM2
  }else{
    map2 <- map2 / NUM2 * NUM1
  }
}
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
RESOLUTION <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1

if(as.character(opt["chromosome"]) != "NA"){
  CHR <- as.character(opt["chromosome"])
}else{
  CHR <- as.character(LocMatrix[1,1])
}
INTRA <- which((as.character(LocMatrix[,1]) == CHR))
LocMatrix <- LocMatrix[INTRA,]

MAX_DISTANCE <- as.numeric(as.character(opt["max_distance"]))
if(MAX_DISTANCE > max(as.numeric(LocMatrix[,3]))){
  MAX_DISTANCE <- max(as.numeric(LocMatrix[,3]))
}

map1 <- cbind(as.character(LocMatrix[INTRA,1]), as.numeric(LocMatrix[INTRA,2]), as.numeric(LocMatrix[INTRA,3])+1, data.frame(map1[INTRA,INTRA]))
map2 <- cbind(as.character(LocMatrix[INTRA,1]), as.numeric(LocMatrix[INTRA,2]), as.numeric(LocMatrix[INTRA,3])+1, data.frame(map2[INTRA,INTRA]))
colnames(map1) <- NULL
colnames(map2) <- NULL


#============================================
# HiC rep program
#============================================
sink("/dev/null")

# check optimal smoothing parameter
h_hat <- htrain(map1, map2, RESOLUTION, MAX_DISTANCE, 0:2)

# pre-processing
Pre_HiC <- prep(map1, map2, RESOLUTION, h_hat, MAX_DISTANCE)

# SCC (Stratum-adjusted Correlation Coefficent) calculation
SCC.out = get.scc(Pre_HiC, RESOLUTION, MAX_DISTANCE)

# SCC score
score <- SCC.out[[3]]

sink()

cat(score, "\n")
