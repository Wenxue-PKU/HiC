#!/usr/bin/Rscript
# Extract target score from Hi-C matrices and output Z-score, distance normalized score

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="target combination file"),
  make_option(c("-m", "--mat"), default="NA", help="matrices"),
  make_option(c("-o", "--out"),help="output files of scores"),
  make_option(c("--format"), default="rds", help="input file format (default: rds)")
)
opt <- parse_args(OptionParser(option_list=option_list))

#=============================================
# test file
#=============================================
# FILE_matrix <- "T:/Project/018_20191128_Noma_Senescent3/data/bmix_backup/IMR90_G_bmix/10kb/ICE2/chr20.rds"
# FILE_in <- "T:/Genome/data/FANTOM5/window/window10kb_combination_chr20.txt"
# FILE_format <- "rds"

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(pbapply)))


options(scipen=10)
FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
FILE_matrix <- as.character(opt["mat"])
FILE_format <- as.character(opt["format"])

if(FILE_format == "rds"){
  FILE_object <- sub(".matrix", ".rds", FILE_matrix)
  if(!file.exists(FILE_object)){
    map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  }else{
    map <- readRDS(FILE_object)
  }
}else if(FILE_format == "matrix"){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  cat("Unknown input file format")
  q()
}
r <- rownames(map)
tmp <- unlist(strsplit(r[1], split = ":"))
RESOLUTION <- as.numeric(tmp[3]) - as.numeric(tmp[2]) + 1

D_table <- fread(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE) %>% as.data.frame()
COLUMN_NAME <- colnames(D_table)
D_table <- D_table[,1:2]
colnames(D_table) <- c("window_left", "window_right")
D_table <- D_table %>% mutate(t1=window_left, t2=window_right) %>% tidyr::separate(t1, c("chr1", "start1", "end1"), ":") %>% tidyr::separate(t2, c("chr2", "start2", "end2"), ":")
D_table <- D_table %>% mutate(distance = abs(as.numeric(start1) - as.numeric(start2))) %>% select(window_left, window_right, distance)


#=============================================
# Calculate average and sd
#=============================================
total_per_line <- apply(map, 1, sum)
NUM_LINE <- nrow(map)
getStat <- function(d){
  distance <- d * RESOLUTION
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  scores <- as.numeric(map[index3])
  scores <- scores[!is.na(scores)]
  if(length(scores) > 2){
    data.frame(distance=distance, ave=mean(scores), sd=sd(scores))
  }else{
    data.frame(distance=distance, ave=NA, sd=NA)
  }
}
MAX_distance <- D_table %>% pull(distance) %>% max()
D_stat <- do.call(rbind, pblapply(1:(MAX_distance/RESOLUTION), getStat))

### filter window
D_table <- D_table %>% filter(window_left %in% r, window_right %in% r)
D_table <- D_table %>% mutate(raw=map[D_table %>% select(window_left, window_right) %>% as.matrix()])
D_table <- dplyr::left_join(D_table, D_stat, by="distance")
rm(D_stat)
D_table <- D_table %>% mutate(disScore=raw/ave, Zscore=(raw-ave)/sd)
D_table <- D_table %>% select(-ave, -sd)
colnames(D_table)[1:2] <- COLUMN_NAME[1:2]


write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


