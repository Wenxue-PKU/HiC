#!/usr/bin/Rscript
# chromosomeとdistanceから対応するaverageのデータを抽出する

### 入力ファイルは、カラム名に少なくとも、chrとdistanceが必要。

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="input file"),
  make_option(c("-d", "--distance"), default="NA", help="distance file"),
  make_option(c("-o", "--out"), default="NA", help="output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))


# ### test data
# FILE_dis <- "X:/hideki_projects/2019-01-18_HiC_human_senescent2/data/IMR90_OIS_mixAll2_distance.txt"
# FILE_in <- "X:/hideki_projects/2019-01-18_HiC_human_senescent2/out/2019-05-01_compare_hic_score_of_EPSP_between_samples/region/OIS_SP_distance.txt"

FILE_in <- as.character(opt["in"])
FILE_dis <- as.character(opt["distance"])
FILE_out <- as.character(opt["out"])


#=============================================
# Distance fileを読み込む
#=============================================
D_dis <- read.table(FILE_dis, header=TRUE, sep="\t", stringsAsFactors = FALSE)
D_dis <- D_dis %>% filter(chromosome1 == chromosome2) %>% select(chromosome1, distance, average)
D_dis <- D_dis %>% rename(chr=chromosome1)


#=============================================
# average取得
#=============================================
D_table <- read.table(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE)
getAverage <- function(i){
  D_sub <- D_dis %>% filter(chr==D_table[i,"chr"]) %>% mutate(len=abs(D_table[i,"distance"]-distance))
  score <- D_sub %>% filter(len==min(D_sub$len)) %>% pull(average) %>% mean()
  score
}
D_out <- data.frame(average=sapply(1:nrow(D_table), getAverage))


write.table(D_out, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



