# BUTLR用に変換する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("--chrom"),help="chrom size file"),
  make_option(c("-o", "--out"),help="output png file")
)
opt <- parse_args(OptionParser(option_list=option_list))

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

# create new column name
chromSize <- read.table(as.character(opt["chrom"]), sep="\t", header=F, check.names = FALSE)
MAX <- chromSize[chromSize[,1]==as.character(LocMatrix[1,1]),2]
Bins <- seq(0, MAX, by=Resolution)
if(max(Bins) + Resolution -1 < MAX){
  Bins <- c(Bins, max(Bins)+Resolution)
}


options(scipen=10)
r.new <- paste(as.character(LocMatrix[1,1]), Bins, Bins+Resolution-1, sep=":")

map.new <- matrix(0, nrow=length(r.new), ncol=length(r.new))
colnames(map.new) <- r.new
rownames(map.new) <- r.new
map.new[r,r] <- map


write.table(map.new, file=as.character(opt["out"]), quote=FALSE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE)


