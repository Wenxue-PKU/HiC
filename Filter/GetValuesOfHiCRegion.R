# 指定したHiC regionsのスコアを出力

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="rds file"),
  make_option(c("-b", "--region"),help="HiC region file"),
  make_option(c("-o", "--out"),help="output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_out <- as.character(opt["out"])
if(file.exists(FILE_out)){
  dummy <- file.remove(FILE_out)
}


# FILE_object <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/data/SCIMR/40kb/Raw/chr1.rds"
FILE_object <- as.character(opt["in"])
map <- readRDS(FILE_object)
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

# FILE_region <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2017-12-02_all_combination_EP_TSS/EP_TSS_region.txt"

FILE_region <- as.character(opt["region"])
con <- file(FILE_region, "r")
while(length(lines <- readLines(con, n=100,warn = F)) > 0){
  REGION <- strsplit(lines, "\t")
  REGION <- matrix(unlist(REGION), ncol=8, byrow=TRUE)
  colnames(REGION) <- c("chr1", "start1", "end1", "id1", "chr2", "start2", "end2", "id2")
  
  getScore <- function(i){
    CHR1 <- as.character(REGION[i,"chr1"])
    START1 <- as.numeric(REGION[i,"start1"])
    END1 <- as.numeric(REGION[i,"end1"])
    CHR2 <- as.character(REGION[i,"chr2"])
    START2 <- as.numeric(REGION[i,"start2"])
    END2 <- as.numeric(REGION[i,"end2"])
    
    index1 <- which((as.character(LocMatrix[,1]) == CHR1) & (as.numeric(LocMatrix[,3]) >= START1) & (as.numeric(LocMatrix[,2]) <= END1))
    index2 <- which((as.character(LocMatrix[,1]) == CHR2) & (as.numeric(LocMatrix[,3]) >= START2) & (as.numeric(LocMatrix[,2]) <= END2))
    max(map[index1, index2], na.rm = T)
  }
  
  Score <- sapply(1:nrow(REGION), getScore)
  OUT <- cbind(REGION, Score)
  write.table(OUT, file = FILE_out, sep="\t", eol = "\n", col.names = F, row.names = F, append = T, quote = F)
}
close(con)






