# HiC score distribution check

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output png file"),
  make_option(c("--max_distance"), default=1000000000000, help="maximum distance to calculate"),
  make_option(c("--title"), default="",help="matrix file")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_matrix <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/data/IMR90_RS2/500kb/Raw/chr1.rds"
FILE_matrix <- "X:/hideki_projects/2017-11-15_HiC_rao_kb/data/GM12878_mix/500kb/Raw/chr1.rds"

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

RESOLUTION <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1

MAX_DISTANCE = as.numeric(opt["max_distance"])
if(MAX_DISTANCE >= length(r) * RESOLUTION){
  MAX_DISTANCE <- (length(r)-1) * RESOLUTION
}

### もしも１列全部が0だったら除く
LINE_SUM <- apply(map, 1, sum, na.rm=TRUE)
ZERO_LINE <- which(LINE_SUM == 0)
map[ZERO_LINE, ZERO_LINE] <- NA


# 指定した距離までのスコアを取得
score <- c()
for(d in 1:as.integer(MAX_DISTANCE / RESOLUTION)){
  index1 <- 1:(nrow(map) - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  score <- c(score, as.numeric(map[index3]))
}
score <- round(score[!is.na(score)], digits = 1)

range <- 0:50
percent <- sapply(range, function(i){sum(score==i)})/length(score)*100


FILE_OUT <- as.character(opt["out"])
png(file=FILE_OUT, width=400, height=300, units="px", bg="white")
par(oma=c(0,0,0,0), mar=c(4,4,3,1))
plot(range, percent, xlab="Read for each combination", ylab="% population", main=as.character(opt["title"]), cex.lab=1.4, cex.axis=1.2, cex.main=1.3)
dummy <- dev.off()





