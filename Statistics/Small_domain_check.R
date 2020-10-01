#!/usr/bin/Rscript
# small domainの存在をチェックするプログラム

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrices file"),
  make_option(c("-o", "--out"),help="output file"),
  make_option(c("--max_distance"), default="100000", help="maximum distance to consider")
)
opt <- parse_args(OptionParser(option_list=option_list))

options(scipen=10)
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(mixtools)))


FILE_in <- as.character(opt["in"])
map <- readRDS(FILE_in)
MAXIMUM_DISTANCE <- as.numeric(as.character(opt["max_distance"]))

map <- ifelse(is.infinite(map), NA, map)
r <- rownames(map)
LocMatrix <- as.data.frame(r) %>% mutate(name=r) %>% tidyr::separate(r, c("chr", "start", "end"), ":", convert=TRUE)

### 複数のchromosomeがある場合、一番長いやつだけにする
if(LocMatrix %>% distinct(chr) %>% nrow() > 1){
  longest_chr <- LocMatrix %>% group_by(chr) %>% tally() %>% arrange(desc(n)) %>% head(n=1) %>% pull(chr)
  LocMatrix <- LocMatrix %>% filter(chr==longest_chr)
  index_chr <- LocMatrix %>% pull(name)
  map <- map[index_chr, index_chr]
}


MAP_RESOLUTION <- LocMatrix %>% head(n=1) %>% mutate(reso=end-start+1) %>% pull(reso)
NUM_LINE <- nrow(map)
MAX_distance <- min(NUM_LINE-1, as.integer(MAXIMUM_DISTANCE / MAP_RESOLUTION))
D_table <- NULL
for(d in 0:MAX_distance){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(LocMatrix[index1, "name"], LocMatrix[index2, "name"])
  score <- as.numeric(map[index3])
  score <- sort(score[!is.na(score)])
  fit <- normalmixEM(score, k = 2) #try to fit two Gaussians
  if(fit$lambda[1] > fit$lambda[2]){
    n1 <- fit$lambda[1]
    n2 <- fit$lambda[2]
  }else{
    n2 <- fit$lambda[1]
    n1 <- fit$lambda[2]
  }
  df <- data.frame(distance=d*MAP_RESOLUTION, score=1-(n1-n2)/n1)
  D_table <- rbind(D_table, df)
}

Average <- D_table %>% filter(distance > 10000, distance < 100000) %>% pull(score) %>% mean()
cat(Average, "\n")
# D_out <- data.frame(input=FILE_in, score=Average)
FILE_out <- as.character(opt["out"])
write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





#=============================================
# 原理を説明するためのグラフ作成
#=============================================
drawTestGraph <- function(){
  FILE_in <- "T:/Project/020_20200207_Noma_HiC_pombe/data/HiC_Double-MHM/5kb/ICE/ALL.rds"
  map <- readRDS(FILE_in)
  map <- ifelse(is.infinite(map), NA, map)
  r <- rownames(map)
  LocMatrix <- as.data.frame(r) %>% mutate(name=r) %>% tidyr::separate(r, c("chr", "start", "end"), ":", convert=TRUE)
  
  ### 複数のchromosomeがある場合、一番長いやつだけにする
  if(LocMatrix %>% distinct(chr) %>% nrow() > 1){
    longest_chr <- LocMatrix %>% group_by(chr) %>% tally() %>% arrange(desc(n)) %>% head(n=1) %>% pull(chr)
    LocMatrix <- LocMatrix %>% filter(chr==longest_chr)
    index_chr <- LocMatrix %>% pull(name)
    map <- map[index_chr, index_chr]
  }
  
  MAP_RESOLUTION <- LocMatrix %>% head(n=1) %>% mutate(reso=end-start+1) %>% pull(reso)
  NUM_LINE <- nrow(map)
  
  CHECK_DISTANCE <- 120000
  d <- as.integer(CHECK_DISTANCE / MAP_RESOLUTION)
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(LocMatrix[index1, "name"], LocMatrix[index2, "name"])
  score <- as.numeric(map[index3])
  score <- sort(score[!is.na(score)])
  fit <- normalmixEM(score, k = 2) #try to fit two Gaussians
  
  ### graph check
  hist(score, freq = FALSE, breaks = 40, main = paste0(CHECK_DISTANCE/1000, " kb distance"), xlab="HiC score")
  lines(score,fit$lambda[1]*dnorm(score,fit$mu[1],fit$sigma[1]), col = "darkgreen")
  lines(score,fit$lambda[2]*dnorm(score,fit$mu[2],fit$sigma[2]), col = "red")
  rug(x=fit$mu, lwd=1.2)
  
  # D_table.wt <- D_table
  # D_table.rad21 <- D_table
  # D_table <- rbind(D_table.wt %>% mutate(sample="wt"), D_table.rad21 %>% mutate(sample="rad21"))
  # suppressWarnings(suppressMessages(library(ggplot2)))
  # suppressWarnings(suppressMessages(library(ggsci)))
  # D_graph <- D_table %>% mutate(diff=abs(log2(lambda1/lambda2)))
  # ggplot(D_graph, aes(distance/1000, diff, col=sample, group=sample))+
  #   geom_line() +
  #   scale_color_simpsons(name="Sample")+
  #   # coord_cartesian(xlim=c(0,500)) +
  #   theme_bw()+
  #   labs(x="Distance (kb)", y="log2(ratio)") +
  #   theme(
  #     text = element_text(size=18),
  #     legend.position = c(.97, .97),
  #     legend.justification = c("right", "top"),
  #     legend.box.just = "right",
  #     legend.margin = margin(6, 6, 6, 6))
}







