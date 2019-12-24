#!/usr/bin/Rscript
# Extract significant associations from 10kb matrices


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output file"),
  make_option(c("--min"), default="80000", help="minimum distance to check(default : 80kb)"),
  make_option(c("--max"), default="2000000", help="maximum distance to check (default : 2Mb)"),
  make_option(c("--control"), default=1, help="minimum required average read for control region"),
  make_option(c("--local"), default=4, help="local enrichment"),
  make_option(c("--fc"), default=4, help="fc threshold to background of same distance"),
  make_option(c("--FDR"), default=0.01, help="threshold of FDR"),
  make_option(c("--background"), default="4", help="threshold for average background score"),
  make_option(c("--all"), default='NULL', help="file name of all scores. (default:NULL = not output)")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))


FILE_out <- as.character(opt["out"])
FILE_all <- as.character(opt["all"])
FILE_in <- as.character(opt["in"])
T_max <- as.numeric(as.character(opt["max"]))
T_min <- as.numeric(as.character(opt["min"]))
T_control <- as.numeric(as.character(opt["control"]))
T_local <- as.numeric(as.character(opt["local"]))
T_fc <- as.numeric(as.character(opt["fc"]))
FDR <- as.numeric(as.character(opt["FDR"]))
T_back <- as.numeric(as.character(opt["background"]))


if(grepl(".rds", FILE_in)){
  map <- readRDS(FILE_in)
}else{
  map <- as.matrix(read.table(FILE_in, header=TRUE, check.names = FALSE))
}


LocList <- strsplit(rownames(map), ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
NUM_LINE <- nrow(map)
Location <- data.frame(id=1:NUM_LINE, chr=LocMatrix[,1], start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]), stringsAsFactors = FALSE)
Resolution <- Location[1,"end"] - Location[1,"start"] + 1
Sliding_window <- Location[2,"start"] - Location[1,"start"]
T_max <- min(c(T_max,  max(Location[,"end"])- Location[1,"start"]))
rm(LocList, LocMatrix)

# remove inter-chromosome
IntraChr <- function(x){
  Location[x,"chr"] == Location[,"chr"]
}
mask_intra <- cbind(sapply(1:NUM_LINE, IntraChr))
map[!mask_intra] <- NA
rm(mask_intra)



#=============================================
# 周りのスコアを見てフィルタリング
#=============================================
mask <- data.frame(
          V1 = c(8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V2 = c(8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V3 = c(8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V4 = c(6, 6, 6, -4, -4, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V5 = c(6, 6, 6, -4, -4, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V6 = c(6, 6, 6, -4, -4, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V7 = c(4, 4, 4, -2, -2, -2, -10, -10, -10, 0, 0, 0, 0, 0, 0),
          V8 = c(4, 4, 4, -2, -2, -2, -10, -20, -10, 0, 0, 0, 0, 0, 0),
          V9 = c(4, 4, 4, -2, -2, -2, -10, -10, -10, 0, 0, 0, 0, 0, 0),
         V10 = c(2, 2, 2, -1, -1, -1, -3, -3, -3, -5, -5, -5, 0, 0, 0),
         V11 = c(2, 2, 2, -1, -1, -1, -3, -3, -3, -5, -5, -5, 0, 0, 0),
         V12 = c(2, 2, 2, -1, -1, -1, -3, -3, -3, -5, -5, -5, 0, 0, 0),
         V13 = c(1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9),
         V14 = c(1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9),
         V15 = c(1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9)
)
mask <- as.matrix(mask)
SHIFT_TO_CENTER <- (nrow(mask)-1)/2 + 1

mask_all <- which(mask > -20, arr.ind = TRUE) - SHIFT_TO_CENTER
mask_center <- which(mask < -9, arr.ind = TRUE) - SHIFT_TO_CENTER
mask_surround <- list()
for(bb in 1:9){
  mask_surround[[as.character(bb)]] <- which(mask==bb, arr.ind = TRUE) - SHIFT_TO_CENTER
}
checkLocalFc <- function(i,j){
  ### 中心をずらす
  index_all <- cbind(mask_all[,1]+i, mask_all[,2]+j)
  index_center <- cbind(mask_center[,1]+i, mask_center[,2]+j)
  sss_all <- as.numeric(map[index_all])
  sss_center <- as.numeric(map[index_center])
  
  ### 中心はNAでない
  if(is.na(map[i,j])){
    return(0)
  }
  
  ### NAは３つ以下
  if(sum(is.na(sss_all)) > 2){
    return(0)
  }
  
  ### 中心が右上の三角形で一番スコアが高い
  if(map[i,j] < max(sss_all, na.rm = TRUE)){
    return(0)
  }
  
  ### 周りの3x3のどのタイルよりも中心3x3のスコアの方がlocal ratio以上高い
  sss_center_ave <- mean(sss_center, na.rm = TRUE)
  if(sss_center_ave == 0){
    return(0)
  }
  fc <- c()
  for(bb in 1:9){
    index_surrounded <- cbind(mask_surround[[as.character(bb)]][,1]+i, mask_surround[[as.character(bb)]][,2]+j)
    average_surrounded <- mean(as.numeric(map[index_surrounded]), na.rm = TRUE)
    if(average_surrounded == 0){
      fold_change <- sss_center_ave
    }else{
      fold_change <- sss_center_ave / average_surrounded
    }
    fc <- c(fc, fold_change)
  }
  ### 最も低いfold-changeを出力
  min(fc)
}


#### 周りのupstream 200kb, downstream 200kbのスコア
checkBackground <- function(i,j){
  x_min <- max(1, i-200000/Resolution)
  x_max <- min(NUM_LINE, i+200000/Resolution)
  y_min <- max(1, j-200000/Resolution)
  y_max <- min(NUM_LINE, j+200000/Resolution)
  score <- as.numeric(map[x_min:x_max, y_min:y_max])
  if(!is.na(map[i,j])){
    back_s <- mean(score, na.rm = TRUE)
  }else{
    back_s <- (sum(score, na.rm = TRUE) - map[i,j])/ (length(score) - 1)
  }
  back_s
}

D_table <- c()
MAX_NUM <- min(c(NUM_LINE-SHIFT_TO_CENTER-1, max(which(Location$end < T_max))-SHIFT_TO_CENTER))
MIN_NUM <- max(c(SHIFT_TO_CENTER + 1, min(which(Location$start > T_min))))
for(d in MIN_NUM:MAX_NUM){
  if((NUM_LINE - d - SHIFT_TO_CENTER) < (SHIFT_TO_CENTER + 1)){
    next
  }
  index1 <- (SHIFT_TO_CENTER + 1):(NUM_LINE - d - SHIFT_TO_CENTER)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  
  ### local fold-change
  wrapper <- function(k){
    checkLocalFc(index3[k,1], index3[k,2])
  }
  local_fc <- sapply(1:nrow(index3), wrapper)
  
  ### background average
  wrapper2 <- function(k){
    checkBackground(index3[k,1], index3[k,2])
  }
  background_score <- sapply(1:nrow(index3), wrapper2)
  
  scores <- as.numeric(map[index3])
  
  ### top1%のデータだけbackgroundから取り除く
  scores_back <- scores[!is.na(scores)]
  scores_back <- scores_back[scores_back < quantile(scores_back, prob=0.99)]
  
  Average <- mean(scores_back)
  if(!is.na(Average) & Average != 0){
    distance <- d*Sliding_window
    Pvalues <- pnorm(scores, mean=Average, sd=sd(scores_back), lower.tail = FALSE)
    df <- data.frame(id1=index1, id2=index2, distance=distance, score=scores, control=Average, 
                     dis_fc=scores/Average, pval=Pvalues, local_fc=local_fc, back_ave=background_score, stringsAsFactors = FALSE)
    df <- df %>% filter(!is.na(Pvalues))
    D_table <- rbind(D_table, df)
  }
}
rm(df, index1, index2, index3, local_fc, background_score, scores, scores_back, distance, Pvalues, map)

D_table <- D_table %>% mutate(qval=p.adjust(pval, method = "BH"))
D_table <- dplyr::left_join(D_table, Location %>% rename(chr1=chr, start1=start, end1=end), by=c("id1"="id"))
D_table <- dplyr::left_join(D_table, Location %>% rename(chr2=chr, start2=start, end2=end), by=c("id2"="id"))
D_table <- D_table[,c("chr1", "start1", "end1", "chr2", "start2", "end2", "distance", "score", "control", "dis_fc", "pval", "qval", "local_fc", "back_ave")]
D_table <- D_table %>% arrange(desc(local_fc))


### output all scores
if(FILE_all != "NULL"){
  write.table(D_table, FILE_all, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


### Filtering
D_table <- D_table %>% filter(dis_fc >= T_fc, control >= T_control, qval <= FDR, local_fc >= T_local, back_ave >= T_back)
write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


