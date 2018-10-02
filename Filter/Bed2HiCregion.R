# bedファイルを入力として組み合わせを出力する
# distanceのフィルターをかける

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="bed file"),
  make_option(c("--min"), default = -1, help="minimum distance for pairs"),
  make_option(c("--max"), default = 500000, help="maximum distance for pairs"),
  make_option(c("-o", "--out"), help="output region definition file")
)
opt <- parse_args(OptionParser(option_list=option_list))



FILE_OUT <- as.character(opt["out"])
FILE_BED <- as.character(opt["in"])
DATA_head <- read.table(FILE_BED, header=F, nrows = 5, stringsAsFactors=F)
classes <- sapply(DATA_head, class)

if(length(classes) > 3){
  classes[4:length(classes)] <- "NULL"
}
BED <- read.table(FILE_BED, header=F, sep="\t", colClasses = classes)
Middle <- as.integer((as.numeric(BED[,2]) + as.numeric(BED[,3]))/2)
BED <- cbind(BED, Middle)
colnames(BED) <- c("chr", "start", "end", "middle")

MIN <- as.numeric(as.character(opt["min"]))
MAX <- as.numeric(as.character(opt["max"]))
NUM <- nrow(BED)

if(file.exists(FILE_OUT)){
  dummy <- file.remove(FILE_OUT)
}

getComb <- function(i){
  TARGET_BED <- BED[i:NUM,]
  distance <- abs(as.numeric(BED[i,"middle"]) - as.numeric(TARGET_BED[,"middle"]))
  index <- which(as.character(BED[i,"chr"]) == as.character(TARGET_BED[,"chr"]) & distance > MIN & distance <= MAX)
  OUT <- cbind(BED[rep(i,length(index)),c("chr", "start", "end")], TARGET_BED[index,c("chr", "start", "end")])
  write.table(OUT, FILE_OUT, sep="\t", quote = F, eol = "\n", row.names = F, col.names = F, append = TRUE)
}

dummy <- sapply(1:NUM, getComb)










