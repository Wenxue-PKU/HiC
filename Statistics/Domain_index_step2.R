#!/usr/bin/Rscript
# Calculate Domain index score


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="Data file for domain index calcualtion. separated with ,"),
  make_option(c("-o", "--out"), default="NA", help="output file prefix"),
  make_option(c("-n", "--names"), default="NA", help="name list of samples"),
  make_option(c("-s", "--size"), default="big", help="'big' domain or 'small' domains")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_ins <- unlist(strsplit(as.character(opt["in"]), ","))
SAMPLE_NAMES <- unlist(strsplit(as.character(opt["names"]), ","))
FILE_prefix <- as.character(opt["out"])
DOMAIN_SIZE <- as.character(opt["size"])

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(ggsci)))
suppressWarnings(suppressMessages(library(scales)))
# suppressWarnings(suppressMessages(library(ggpubr)))
# suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(pbapply)))


if(length(FILE_ins) != length(SAMPLE_NAMES)){
  cat("Sample names and input file length should be same")
  q()
}


if(DOMAIN_SIZE == "big"){
  Threshold <- 0.9
}else{
  Threshold <- 0
}


readCondition <- function(i){
  n = SAMPLE_NAMES[i]
  file <- FILE_ins[i]
  df <- fread(file, header=TRUE)
  df1 <- df %>% filter(doNR==doNL, doNR != -1)
  df1 <- df1 %>% group_by(doNR) %>% summarize(score1=mean(NormScore)) %>% as.data.frame()
  df2 <- df %>% filter(doNR-doNL==1, doNR != -1)
  df2 <- df2 %>% group_by(doNR) %>% summarize(score2=mean(NormScore)) %>% as.data.frame()
  df3 <- df %>% filter(doNR-doNL==1, doNR != -1)
  df3 <- df3 %>% group_by(doNL) %>% summarize(score3=mean(NormScore)) %>% as.data.frame()
  df <- dplyr::left_join(df1, df2, by=c("doNR"))
  df <- dplyr::left_join(df, df3, by=c("doNR"="doNL"))
  df <- df %>% mutate(control=apply(df[,c("score2", "score3")], 1, mean, na.rm=TRUE))
  df <- df %>% mutate(score=score1/control) 
  df <- df %>% mutate(ranking=percent_rank(score))
  df <- df %>% rename(id=doNR)
  df <- df %>% mutate(name=n) %>% select(name, id, score, ranking)
  df
}
table <- do.call(rbind, pblapply(1:(length(SAMPLE_NAMES)), readCondition))


outputGraph <- function(TITLE, D_table, threshold){
  D_table <- D_table %>% mutate(name=factor(name, levels = SAMPLE_NAMES))
  D_table <- D_table %>% filter(ranking > threshold)
  p <- ggplot(D_table, aes(x=name, y=score)) +
    theme_bw()+
    geom_jitter(alpha = 0.8, position=position_jitter(), size=1.4, col='gray') + 
    geom_boxplot(alpha = 0.3, col="black", outlier.shape = NA) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none", text = element_text(size=14)) +
    stat_summary(fun.y=mean, geom="point", color="green", size=3, alpha=0.6) + 
    labs(y=TITLE, x="")
  save_plot(paste0(FILE_prefix, "_graph.png"), p, base_width = 12, base_height = 8)
  
  ### summary
  D_out <- D_table %>% group_by(name) %>% summarize(score=mean(score)) %>% as.data.frame()
  write.table(D_out, paste0(FILE_prefix, "_index.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
outputGraph("Ave. disNorm. score intra/inter-domains", table, Threshold)






