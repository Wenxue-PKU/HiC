# 指定したbinの長さのsectionを作成する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("--organism"), default="pombe", help="pombe or human or mouse"),
  make_option(c("--resolution"), default="100000", help="bin size"),
  make_option(c("-o", "--out"), help="outpu file name")
)
opt <- parse_args(OptionParser(option_list=option_list))


ORGANISM <- as.character(opt["organism"])
switch(ORGANISM,
       "pombe"={
         CHRs=c("I", "II", "III")
         LENGTH=c(5579133, 4539804, 2452883)
       },
       "human"={
         CHRs=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
         LENGTH=c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
       },
       "mouse"={
         CHRs=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
         LENGTH=c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566, 171031299, 91744698)
       }
)

BIN_SIZE = as.numeric(as.character(opt["resolution"]))

options(scipen=10)
getBins <- function(i){
  Start <- seq(0, LENGTH[i], by=BIN_SIZE)
  End <- Start + BIN_SIZE - 1
  cbind(CHRs[i], Start, End)
}

DATA <- sapply(1:length(CHRs), getBins)

OUT <- c()
for(i in seq(1, length(CHRs))){
  OUT <- rbind(OUT, DATA[[i]])
}

write.table(OUT, file = as.character(opt["out"]), row.names = FALSE, col.names = FALSE, sep="\t", eol = "\n", quote = FALSE)



