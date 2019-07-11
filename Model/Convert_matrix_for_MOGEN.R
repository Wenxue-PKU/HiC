#!/usr/bin/Rscript
# Convert matrix for MOGEN script

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output interaction data file name"),
  make_option(c("--fill"), default="NA", help="fill blank score to this ratio value (ex. 0.5 means top 50% score) "),
  make_option(c("--break_point"), default="NA", help="file name which output breakpoint of chromosomes")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))

FILE_out <- as.character(opt["out"])
FILE_break <- as.character(opt["break_point"])
FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}

LocList <- strsplit(rownames(map), ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
NUM_LINE <- nrow(map)
Location <- data.frame(id=1:NUM_LINE, chr=LocMatrix[,1], start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]), stringsAsFactors = FALSE)
rm(LocList, LocMatrix)

### breakpointを記したファイル
df_break <- Location %>% group_by(chr) %>% summarise(id=max(id)) %>% as.data.frame() %>% mutate(chr_ID=row_number()) %>% select(chr_ID, id)
write.table(df_break, FILE_break, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


colnames(map) <- 1:NUM_LINE
rownames(map) <- 1:NUM_LINE
map <- as.data.frame(map, stringsAsFactors=FALSE)
map <- cbind(1:NUM_LINE, map)
colnames(map)[1] <- "origin"
map <- map %>% tidyr::gather(key="target", value="score", -origin)
map$target <- as.integer(map$target)
map <- map %>% filter(origin < target)

NUM_FILL <- as.character(opt["fill"])
if(NUM_FILL == "NA"){
  map <- map %>% filter(score != 0)
}else{
  T <- quantile(map$score,prob=1-as.numeric(NUM_FILL))
  map <- map %>% mutate(ifelse(score == 0, T, score))
}
map <- map %>% arrange(origin, target)

write.table(map, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)







