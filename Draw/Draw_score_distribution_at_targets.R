# contact mapを描画する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output png file"),
  make_option(c("--title"), default="", help="title of png"),
  make_option(c("--groupname"), default="NA", help="name of group separated by ,"),
  make_option(c("--bed"),help="target bed file")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  saveRDS(map, FILE_object)
}else{
  map <- readRDS(FILE_object)
}

# log2(Observed / Expect)のmatrixに変換する
map_expect <- map
NUM_LINE <- nrow(map)
for(d in 0:(NUM_LINE-1)){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  index4 <- cbind(index2, index1)
  Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
  map_expect[index3] <- Average
  map_expect[index4] <- Average
}
map <- ifelse(map_expect == 0, 0, log2(map / map_expect))
map <- ifelse(is.infinite(map), 0, map)


# read bed file
bed <- as.matrix(read.table(as.character(opt["bed"]), header=FALSE, check.names=FALSE, sep="\t"))
options(scipen=10)
rownames(bed) <- paste(bed[,1], bed[,2], bed[,3], sep=":")
colnames(bed) <- c("chr", "start", "end", "score", "group")

r.common <- intersect(rownames(map), rownames(bed))
# cat("Common section:", length(r.common), "\n")

map <- map[r.common, r.common]
bed <- bed[r.common, ]

LocList <- strsplit(r.common, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)



IntraBed <- function(x){
  score <- rep(NA, length(r.common))
  targetGroup <- bed[x,"group"]
  con <- ( targetGroup == bed[,"group"] & (as.character(bed[x,"chr"]) == as.character(bed[, "chr"])) )
  score[con] <- targetGroup
  score
}

mask <- cbind(sapply(1:length(r.common), IntraBed))
rownames(mask) <- r.common
colnames(mask) <- r.common


Groups <- sort(unique(bed[,"group"]))
if(as.character(opt["groupname"]) != "NA"){
  Groupnames <- unlist(strsplit(as.character(opt["groupname"]), ","))
}else{
  Groupnames <- Groups
}
DATA <- list()
for(g in Groups){
  DATA[[g]] <- as.numeric(map[mask==g])
}

if( length(Groups) == 2){

  xmin <- min(map, na.rm=TRUE)
  xmax <- max(map, na.rm=TRUE)
  
  # cat("xmin", xmin, "\n")
  # cat("xmax", xmax, "\n")
  
  b <- seq(xmin-1, xmax+1, by=0.1)
  
  png(filename=as.character(opt["out"]), width=500, height=400,  units="px", bg="white")
  par(oma=c(0,0,0,0), mar=c(5,5,3,1))
  hist(DATA[[Groups[1]]], col=rgb(0.9,0,0,0.5), xlim=c(-3,3), ylim=c(0,0.5), breaks = b, main=as.character(opt["title"]), 
       xlab="log2(Obs/Exp)", probability = T, cex.axis=1.8, cex.lab=2, cex.main=1.8)
  hist(DATA[[Groups[2]]], col=rgb(0,0,0.9,0.5), add=T, breaks = b, probability = T, cex.axis=1.8, cex.lab=2, cex.main=1.8)
  legend("topright", legend=Groupnames, col= c(rgb(0.9,0,0,0.5), rgb(0,0,0.9,0.5)), 
         lty=1, lwd=3, bty='n', cex=1.5)
  dummy <- dev.off()
}








