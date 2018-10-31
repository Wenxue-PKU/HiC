#!/usr/bin/Rscript
# make combination of two bed files (filter by distance)

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="bed file"),
  make_option(c("--min"), default = -1, help="minimum distance for pairs"),
  make_option(c("--max"), default = 500000, help="maximum distance for pairs"),
  make_option(c("-o", "--out"), help="output region definition file")
)
opt <- parse_args(OptionParser(option_list=option_list))

### Input BED file format
# 1: chr
# 2: start
# 3: end
# 4: region name (option)
# 5-: ignored


FILE_OUT <- as.character(opt["out"])
FILE_BED <- as.character(opt["in"])
DATA_head <- read.table(FILE_BED, header=F, nrows = 5, stringsAsFactors=F)
classes <- sapply(DATA_head, class)

if(length(classes) > 4){
  classes[5:length(classes)] <- "NULL"
}
BED <- read.table(FILE_BED, header=F, sep="\t", colClasses = classes)
if(length(classes) > 3){
  row_index <- BED[,4]
}else{
  row_index <- 1:nrow(BED)
}
Middle <- as.integer((as.numeric(BED[,2]) + as.numeric(BED[,3]))/2)
BED <- cbind(BED[,1], BED[,2], BED[,3], Middle, row_index)
colnames(BED) <- c("chr", "start", "end", "middle", "row_index")

MIN <- as.numeric(as.character(opt["min"]))
MAX <- as.numeric(as.character(opt["max"]))
NUM <- nrow(BED)

if(file.exists(FILE_OUT)){
  dummy <- file.remove(FILE_OUT)
}

# Output is always left side is smaller
getComb <- function(i){
  TARGET_BED <- BED[i:NUM,]
  distance <- abs(as.numeric(BED[i,"middle"]) - as.numeric(TARGET_BED[,"middle"]))
  index <- which(as.character(BED[i,"chr"]) == as.character(TARGET_BED[,"chr"]) & distance > MIN & distance <= MAX)
  OUT <- cbind(BED[rep(i,length(index)),c("chr", "start", "end", "row_index")], TARGET_BED[index,c("chr", "start", "end", "row_index")])
  index <- which(as.numeric(OUT[,2]) > as.numeric(OUT[,6])) 
  OUT[index,] <- OUT[index,c(5,6,7,8,1,2,3,4)]
  write.table(OUT, FILE_OUT, sep="\t", quote = F, eol = "\n", row.names = F, col.names = F, append = TRUE)
}

dummy <- sapply(1:NUM, getComb)










