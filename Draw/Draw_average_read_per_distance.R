#!/usr/bin/Rscript
# Distanceごとのreadの分布のグラフを作成する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output png prefix"),
  make_option(c("--min"), default="20000", help="minimum distance to check(default : 20kb)"),
  make_option(c("--max"), default="1000000", help="maximum distance to check (default : 1Mb)"),
  make_option(c("--title"), default="", help="title of graph")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))

OUT_PREFIX <- as.character(opt["out"])
FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}

LocList <- strsplit(rownames(map), ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
NUM_LINE <- nrow(map)
T_max <- as.numeric(as.character(opt["max"]))
T_min <- as.numeric(as.character(opt["min"]))
TITLE <- as.character(opt["title"])

# inter-chromosomeは除外する
IntraChr <- function(x){
  as.character(LocMatrix[x,1]) == as.character(LocMatrix[,1])
}
mask_intra <- cbind(sapply(1:NUM_LINE, IntraChr))
map[!mask_intra] <- NA
rm(mask_intra, LocMatrix, LocList)

D_table <- c()
D_qvals <- c()
MIN_sigs <- c()
NUM_sigs <- c()
NUM_comb <- c()
MAX_NUM <- min(c(NUM_LINE-1, as.integer(T_max/Resolution)))
MIN_NUM <- max(c(1, as.integer(T_min/Resolution)))
target <- seq(MIN_NUM, MAX_NUM, by=5)
for(d in MIN_NUM:MAX_NUM){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  
  scores <- as.numeric(map[index3])
  scores <- scores[!is.na(scores)]
  
  Average <- mean(scores)
  if(!is.na(Average) & Average != 0){
    bw <- max(1, (max(scores) - min(scores))/30)
    dn <- density(scores, bw=bw)
    distance <- d*Resolution/1000
    Pvalues <- pnorm(scores, mean=mean(scores), sd=sd(scores), lower.tail = FALSE)
    qval <- p.adjust(Pvalues, method = "BH")
    NUM_sigs <- c(NUM_sigs, sum(qval < 0.05))
    NUM_comb <- c(NUM_comb, length(scores))
    
    if(d %in% target){
      df <- data.frame(distance=distance, read=dn$x, density=dn$y, stringsAsFactors = FALSE)
      D_table <- rbind(D_table, df)
      MIN_sigs <- c(MIN_sigs, as.integer(min(scores[qval < 0.05], na.rm = TRUE)))
      df <- data.frame(distance=distance, qval=qval, stringsAsFactors = FALSE)
      D_qvals <- rbind(D_qvals, df)
    }
  }
}
rm(df, dn, index3, distance, Pvalues, qval)

c <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                        "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(D_table %>% distinct(distance) %>% nrow())


D_table <- D_table %>% mutate(distance=factor(distance)) %>% arrange((distance))
xmax <- D_table %>% filter(density > 0.0001) %>% pull(read) %>% max()
ymax <- D_table %>% pull(density) %>% max()
p <- ggplot(D_table, aes(x=read, y=density, col=distance)) + geom_line() + theme_bw() + 
  scale_color_manual(values=c, name="Distance") +
  theme(legend.position="right") +
  coord_cartesian(xlim=c(0, xmax)) +
  geom_rug(data=data.frame(MIN_sigs), aes(x=MIN_sigs), sides="t", inherit.aes = F, col=c, lwd=1.2) +
  annotate("text", x=MIN_sigs, y=ymax, label=MIN_sigs) +
  labs(title=TITLE)
save_plot(paste0(OUT_PREFIX, "_read_distribution_per_distance.png"), p, base_width = 6, base_height = 5)

D_number <- data.frame(distance=MIN_NUM:MAX_NUM*Resolution/1000, sig=NUM_sigs, combination=NUM_comb)
D_number <- D_number %>% mutate(per=sig/combination*100)


p <- ggplot(D_number, aes(x=distance, y=per)) + geom_bar(stat = "identity") + 
  labs(x="Distance(kb)", y="% of significant pairs (FDR < 0.05)", title=TITLE)
save_plot(paste0(OUT_PREFIX, "_percent_of_significant_pair.png"), p, base_width = 5, base_height = 4)




p <- ggplot(D_qvals, aes(x=qval)) + geom_histogram(aes(y=..density..),breaks=seq(0,1, by=0.05), color='black', fill='grey')+
  coord_cartesian(ylim=c(0, 0.5)) +
  facet_wrap(~distance)+labs(x="FDR", Y="Density") +
  labs(title=TITLE)
save_plot(paste0(OUT_PREFIX, "_FDR_distribution.png"), p, base_width = 12, base_height = 10)

