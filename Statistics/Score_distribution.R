# HiC score distribution check

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output text file"),
  make_option(c("--max_distance"), default=1000000000000, help="maximum distance to calculate")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(ggsci)))



FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])

if(grepl(".rds", FILE_in)){
  map <- readRDS(FILE_in)
}else{
  map <- as.matrix(read.table(FILE_in, header=TRUE, check.names = FALSE))
}

r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)

RESOLUTION <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1

MAX_DISTANCE = as.numeric(opt["max_distance"])
if(MAX_DISTANCE >= length(r) * RESOLUTION){
  MAX_DISTANCE <- (length(r)-1) * RESOLUTION
}

### NAは0に置換する
map <- ifelse(is.na(map), 0, map)

# 指定した距離までのスコアを取得
getScore <- function(d){
  index1 <- 1:(nrow(map) - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  data.frame(score=round(map[index3]))
}

D_data <- do.call(rbind, lapply(1:as.integer(MAX_DISTANCE / RESOLUTION), getScore))

D_table <- rbind(D_data %>% filter(score != 0) %>% mutate(cate_sameNum=ntile(score, 19) + 1),
                 D_data %>% filter(score == 0) %>% mutate(cate_sameNum=1))
D_table <- D_table %>% mutate(cate_sameWidth=cut_interval(score, n = 20))


D_summay <- cbind(
  D_table %>% group_by(cate_sameNum) %>% summarize(min=min(score), max=max(score), n=n(), p=format(n()/nrow(D_table) * 100, digits = 3)) %>%
    as.data.frame() %>% mutate(out=paste0(min, "-", max, ": ", n, " (", p, "%)")) %>% select(out),
  
  D_table %>% group_by(cate_sameWidth) %>% summarize(min=min(score), max=max(score), n=n(), p=format(n()/nrow(D_table) * 100, digits = 3)) %>%
    as.data.frame() %>% mutate(out=paste0(min, "-", max, ": ", n, " (", p, "%)")) %>% select(out)
)
colnames(D_summay) <- c("same number", "same width")

write.table(D_summary, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





