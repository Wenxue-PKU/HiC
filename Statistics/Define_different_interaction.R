#!/usr/bin/Rscript
# Define different interactions

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-n", "--names"), default="NA", help="sample names separated by ,"),
  make_option(c("-i", "--in"), default="NA", help="list of matrices separated by ,"),
  make_option(c("-o", "--out"), default="NA", help="output directory. XXX.png and xxx.txt will be output"),
  make_option(c("-w", "--wt"), default="NA", help="two index showing which samples are biological replicates ex. 1,2"),
  make_option(c("-c", "--comparison"), default="NA", help="index showing which combination should be check. Separated with :. first index devided by second index.  Multiple comparison could be separated by ,. ex. 1:2, 1:3"),
  make_option(c("--category_number"), default="50", help="number of grouping for average calcualtion")
)
opt <- parse_args(OptionParser(option_list=option_list))

MAP_LIST <- unlist(strsplit(as.character(opt["in"]), ","))
DIR_out <- as.character(opt["out"])
SAMPLES <- unlist(strsplit(as.character(opt["names"]), ","))
INDEX_WT <- as.integer(unlist(strsplit(as.character(opt["wt"]), ",")))
INDEX_COMB <- unlist(strsplit(as.character(opt["comparison"]), ","))
NUMBER_GROUP <- as.numeric(as.character(opt["category_number"]))

if(length(SAMPLES) != length(MAP_LIST)){
  cat(paste0("map numbers(", length(MAP_LIST), ") and sample number(", length(SAMPLES), ") should be equals\n"))
  q()
}

options(scipen=10)
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(ggsci)))
# devtools::install_github("thomasp85/patchwork")
suppressWarnings(suppressMessages(library(patchwork)))


D_mat <- NULL
for(ss in 1:length(SAMPLES)){
  map <- readRDS(MAP_LIST[ss])
  r <- rownames(map)
  df <- as.data.frame(map) %>% mutate(bin1=r) %>% tidyr::gather(key=bin2, value=score, -bin1)
  df <- df %>% rename(!!as.name(SAMPLES[ss]) := score)
  
  if(is.null(D_mat)){
    D_mat <- df
  }else{
    D_mat <- dplyr::left_join(D_mat, df, by=c("bin1", "bin2"))
  }
}
rm(df, map)


splitLocation <- function(r){
  LocList <- strsplit(r, ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  Location <- data.frame(chr=LocMatrix[,1], start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]), stringsAsFactors = FALSE)
  Location
}
Loc_L <- splitLocation(D_mat$bin1)
Loc_R <- splitLocation(D_mat$bin2)
D_mat <- D_mat %>% mutate(start1=Loc_L$start, end1=Loc_L$end, start2=Loc_R$start, end2=Loc_R$end)
rm(Loc_L, Loc_R)


D_mat <- D_mat %>% mutate(distance=abs(start1-start2)/1000)

### 0を除く
tmp <- apply(D_mat[,SAMPLES], 1, min)
D_mat <- D_mat %>% filter(tmp!=0)
rm(tmp)


D_mat<- D_mat %>% filter(start1 < start2)
D_ratio_replica <- D_mat %>% mutate(ratio=log2(!!as.name(SAMPLES[INDEX_WT[1]])/!!as.name(SAMPLES[INDEX_WT[2]])), 
                                    ave=(!!as.name(SAMPLES[INDEX_WT[1]])+!!as.name(SAMPLES[INDEX_WT[2]]))/2) %>% select(ratio, ave, distance)

### calculate c2/c1
DEanalysis <- function(c2, c1){
  D_ratio <- D_mat %>% mutate(ratio=log2(!!as.name(c2)/!!as.name(c1)), ave=(!!as.name(c1)+!!as.name(c2))/2) %>% select(bin1, bin2, ratio, ave, distance)

  p1 <- ggplot(data = D_ratio) +
    aes(x = distance, y = ratio) +
    geom_point(alpha=0.1) +
    geom_smooth(span = 0.83) +
    scale_color_startrek(name="") +
    labs(x = "Distance between pair (kb)",
         y = paste0("log2 (", c2, " / ", c1, ")"),
         title = "Before normalize") +
    theme_bw() +
    theme(legend.justification=c(1,0), legend.position=c(1,0), text = element_text(size=11))

  
  ### Average of log2 ratio
  D_ratio <- D_ratio %>% mutate(dcate=dplyr::ntile(distance, 50))
  T_factor <- D_ratio %>% group_by(dcate) %>% summarize(ff=mean(ratio)) %>% as.data.frame()
  
  ### replicate ratio
  D_ratio_replica <- D_ratio_replica %>% mutate(dcate=dplyr::ntile(distance, 50)) 
  T_factor_replica <- D_ratio_replica %>% group_by(dcate)  %>% summarize(ff=mean(ratio)) %>% as.data.frame()
  
  
  ### T_factorでnormalize
  D_ratio <- dplyr::left_join(D_ratio, T_factor, by=c("dcate")) %>% mutate(ratio.norm=ratio - ff) %>% select(-ff, -dcate)
  D_ratio_replica <- dplyr::left_join(D_ratio_replica, T_factor_replica, by=c("dcate")) %>% mutate(ratio.norm=ratio - ff) %>% select(-ff, -dcate)
  rm(T_factor, T_factor_replica)
  
  ### normalize後
  p2 <- ggplot(data = D_ratio) +
    aes(x = distance, y = ratio.norm) +
    geom_point(alpha=0.1) +
    geom_smooth(span = 0.83) +
    scale_color_startrek(name="") +
    labs(x = "Distance between pair (kb)",
         y = paste0("log2 (", c2, " / ", c1, ")"),
         title = "After normalize") +
    theme_bw() +
    theme(legend.justification=c(1,0), legend.position=c(1,0), text = element_text(size=11))
  
  
  ### averageをNUMBER_GROUPに区切る
  D_ratio_replica2 <- rbind(D_ratio_replica, D_ratio_replica %>% mutate(ratio.norm=ratio.norm*-1))
  
  D_background <- D_ratio_replica2 %>% mutate(category=dplyr::ntile(ave, NUMBER_GROUP))
  D_background.stat <- D_background %>% group_by(category) %>% summarize(n=n(), min=min(ave), max=max(ave)) %>% as.data.frame()
  
  
  p3 <- ggplot(data = D_ratio_replica2) +
    aes(x = ave, y = ratio.norm) +
    geom_point(color='grey20', alpha=0.1) +
    scale_x_log10() +
    geom_hline(yintercept = 0, col='red') +
    geom_vline(xintercept = D_background.stat %>% pull(max), col='blue') + 
    labs(x = paste0("Average read between ", SAMPLES[INDEX_WT[1]], " and ", SAMPLES[INDEX_WT[2]]),
         y = paste0("log2 (", SAMPLES[INDEX_WT[1]], " / ", SAMPLES[INDEX_WT[2]], ")"),
         title = paste0("Fluctuation of replicate (Grouped ", NUMBER_GROUP, " categories)")) +
    theme_bw() +
    theme(text = element_text(size=11))
  
  
  D_ratio2 <- NULL
  for(cc in D_background.stat %>% pull(category)){
    scores <- D_background %>% filter(category == cc) %>% pull(ratio.norm)
    MIN <- D_background.stat %>% filter(category == cc) %>% pull(min)
    MAX <- D_background.stat %>% filter(category == cc) %>% pull(max)
    
    if(cc == NUMBER_GROUP){
      tmp <- D_ratio %>% filter(ave >= MIN)
    }else{
      tmp <- D_ratio %>% filter(ave >= MIN & ave <= MAX)
    }
    p <- pnorm(tmp[,"ratio.norm"], mean=0, sd=sd(scores), lower.tail = FALSE)
    p <- ifelse(p > 0.5, 1-p, p)
    D_ratio2 <- rbind(D_ratio2, data.frame(tmp, pval=p))
  }
  D_ratio2 <- D_ratio2 %>% mutate(FDR=p.adjust(pval, method = "BH"))
  
  
  ### P-value distribution check
  # hist(D_ratio2 %>% pull(pval))
  result <- D_ratio2 %>% summarize(up=sum(ratio.norm > 0 & FDR < 0.05), down=sum(ratio.norm < 0 & FDR < 0.05))
  
  
  p4 <- ggplot(data = D_ratio2) +
    aes(x = ave, y = ratio.norm) +
    geom_point(color = 'grey20', alpha=0.3) +
    geom_hline(yintercept = 0, col='blue') +
    geom_point(data=subset(D_ratio2, FDR < 0.05), col='red', alpha=0.3) +
    scale_x_log10() +
    labs(x = paste0("Average read between ", c1, " and ", c2),
         y = paste0("log2 (", c2, " / ", c1, ")"),
         title = paste0("Siginificant change ", c2, " / ", c1),
         subtitle = paste0("Up: ", result$up, ", Down: ", result$down)) +
    theme_bw() +
    theme(text = element_text(size=11))
  
  pmix <- (p1 | p2 | p3)/p4
  save_plot(paste0(DIR_out, "Comparison_of_", c1, "_and_", c2, ".png"), pmix, base_height = 10, base_width = 15, dpi=100)

  D_out <- D_ratio2 %>% filter(FDR < 0.05) %>% select(bin1, bin2, ratio.norm, pval, FDR)
  D_out <- D_out %>% rename(log2_fc=ratio.norm)
  D_out <- dplyr::left_join(D_out, D_mat %>% select(bin1, bin2, !!as.name(c1), !!as.name(c2)), by=c("bin1", "bin2"))
  write.table(D_out, paste0(DIR_out, "Comparison_of_", c1, "_and_", c2, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

for(mm in INDEX_COMB){
  cc <- as.integer(unlist(strsplit(mm, ":")))
  DEanalysis(SAMPLES[cc[1]], SAMPLES[cc[2]])
}


