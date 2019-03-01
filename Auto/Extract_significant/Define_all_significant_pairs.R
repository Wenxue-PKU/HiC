#!/usr/bin/Rscript
# Extract significant associations


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output file"),
  make_option(c("--min"), default="20000", help="minimum distance to check(default : 20kb)"),
  make_option(c("--max"), default="2000000", help="maximum distance to check (default : 2Mb)")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))

FILE_out <- as.character(opt["out"])
FILE_matrix <- as.character(opt["in"])
T_max <- as.numeric(as.character(opt["max"]))
T_min <- as.numeric(as.character(opt["min"]))
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
Resolution <- Location[1,"end"] - Location[1,"start"] + 1
rm(LocList, LocMatrix)

# remove inter-chromosome
IntraChr <- function(x){
  Location[x,"chr"] == Location[,"chr"]
}
mask_intra <- cbind(sapply(1:NUM_LINE, IntraChr))
map[!mask_intra] <- NA
rm(mask_intra)


D_table <- c()
MAX_NUM <- min(c(NUM_LINE-1, as.integer(T_max/Resolution)))
MIN_NUM <- max(c(1, as.integer(T_min/Resolution)))
for(d in MIN_NUM:MAX_NUM){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  scores <- as.numeric(map[index3])
  Average <- mean(scores, na.rm=TRUE)
  if(!is.na(Average) & Average != 0){
    distance <- d*Resolution
    Pvalues <- pnorm(scores, mean=Average, sd=sd(scores, na.rm = TRUE), lower.tail = FALSE)
    df <- data.frame(id1=index1, id2=index2, distance=distance, score=scores, control=Average, 
                     fc=log2(scores/Average), pval=Pvalues, stringsAsFactors = FALSE)
    df <- df %>% filter(!is.na(Pvalues))
    D_table <- rbind(D_table, df)
  }
}
rm(df, map, index3)

NUM_all <- D_table %>% nrow()
D_table <- D_table %>% mutate(qval=p.adjust(pval, method = "BH")) %>% filter(qval < 0.05)
NUM_sig <- D_table %>% nrow()
cat(paste(NUM_sig, NUM_all, sep="\t"))

D_table <- dplyr::left_join(D_table, Location %>% rename(chr1=chr, start1=start, end1=end), by=c("id1"="id"))
D_table <- dplyr::left_join(D_table, Location %>% rename(chr2=chr, start2=start, end2=end), by=c("id2"="id"))
D_table <- D_table[,c("chr1", "start1", "end1", "chr2", "start2", "end2", "distance", "score", "control", "fc", "pval", "qval")]
D_table <- D_table %>% arrange(qval)
  
write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






