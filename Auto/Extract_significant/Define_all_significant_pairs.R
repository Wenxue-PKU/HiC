#!/usr/bin/Rscript
# Extract significant associations


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output file"),
  make_option(c("--min"), default="20000", help="minimum distance to check(default : 20kb)"),
  make_option(c("--max"), default="2000000", help="maximum distance to check (default : 2Mb)"),
  make_option(c("--control"), default=1, help="minimum required average read for control region"),
  make_option(c("--FDR"), default=0.01, help="threshold of FDR")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))

FILE_out <- as.character(opt["out"])
FILE_matrix <- as.character(opt["in"])
T_max <- as.numeric(as.character(opt["max"]))
T_min <- as.numeric(as.character(opt["min"]))
T_control <- as.numeric(as.character(opt["control"]))
FDR <- as.numeric(as.character(opt["FDR"]))


FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}

LocList <- strsplit(rownames(map), ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
NUM_LINE <- nrow(map)
Location <- data.frame(id=1:NUM_LINE, chr=LocMatrix[,1], start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]), stringsAsFactors = FALSE)
Resolution <- Location[1,"end"] - Location[1,"start"] + 1
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
          V1 = c(8, 8, 8, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V2 = c(8, 8, 8, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V3 = c(8, 8, 8, 10, 10, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V4 = c(6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V5 = c(6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V6 = c(6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          V7 = c(4, 4, 4, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0),
          V8 = c(4, 4, 4, 0, 0, 0, -1, -2, -1, 0, 0, 0, 0, 0, 0),
          V9 = c(4, 4, 4, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0),
         V10 = c(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 11, 11),
         V11 = c(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 11, 11),
         V12 = c(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 11, 11),
         V13 = c(1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9),
         V14 = c(1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9),
         V15 = c(1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9)
)

SEARCH_AREA <- 7

mask <- as.matrix(mask)
mask_all <- which(mask>-2, arr.ind = TRUE)
mask_center <- which(mask<0, arr.ind = TRUE)
mask_surround <- list()
for(bb in 1:9){
  mask_surround[[as.character(bb)]] <- which(mask==bb, arr.ind = TRUE)
}
checkSurrounded <- function(i,j){
  ### 中心をずらす
  index_all <- cbind(mask_all[,1]+i-SEARCH_AREA -1, mask_all[,2]+j-SEARCH_AREA-1)
  index_center <- cbind(mask_center[,1]+i-SEARCH_AREA-1, mask_center[,2]+j-SEARCH_AREA-1)
  sss_all <- as.numeric(map[index_all])
  sss_center <- as.numeric(map[index_center]) 
  
  ### 中心はNAでない
  if(is.na(map[i,j])){
    return(0)
  }
  
  ### 90% 以上がスコアで埋まっている
  if((sum(!is.na(sss_all)) / length(sss_all)) < 0.90){
    return(0)
  }
  
  ### centerは少なくとも3つ以上欠けていてはダメ
  if(sum(is.na(sss_center)) > 3){
    return(0)
  }
  
  ### 周りの3x3のどのタイルよりも中心3x3のスコアの方が高い
  sss_center_ave <- mean(sss_center, na.rm = TRUE)
  for(bb in 1:9){
    index_surrounded <- cbind(mask_surround[[as.character(bb)]][,1]+i-SEARCH_AREA-1, mask_surround[[as.character(bb)]][,2]+j-SEARCH_AREA-1)
    if(mean(as.numeric(map[index_surrounded]), na.rm = TRUE) > sss_center_ave){
      return(0)
    }
  }
  return(1)
}


D_table <- c()
MAX_NUM <- min(c(NUM_LINE-1-SEARCH_AREA, as.integer(T_max/Resolution)))
MIN_NUM <- max(c(SEARCH_AREA + 1, as.integer(T_min/Resolution)))
for(d in MIN_NUM:MAX_NUM){
  index1 <- (SEARCH_AREA + 1):(NUM_LINE - d - SEARCH_AREA)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  
  wrapper <- function(k){
    checkSurrounded(index3[k,1], index3[k,2])
  }
  Okay <- sapply(1:nrow(index3), wrapper)
  
  scores <- as.numeric(map[index3])
  
  ### top1%のデータだけbackgroundから取り除く
  scores_back <- scores[!is.na(scores)]
  scores_back <- scores_back[scores_back < quantile(scores_back, prob=0.99)]
  
  Average <- mean(scores_back)
  if(!is.na(Average) & Average != 0){
    distance <- d*Resolution
    Pvalues <- pnorm(scores, mean=Average, sd=sd(scores_back), lower.tail = FALSE)
    df <- data.frame(id1=index1, id2=index2, distance=distance, score=scores, control=Average, 
                     fc=log2(scores/Average), pval=Pvalues, okay=Okay, stringsAsFactors = FALSE)
    df <- df %>% filter(!is.na(Pvalues))
    D_table <- rbind(D_table, df)
  }
}
rm(df, index1, index2, index3, Okay, scores, Pvalues)

cat(paste("Total target: ", D_table %>% nrow(), "\n"))
D_table <- D_table %>% mutate(qval=p.adjust(pval, method = "BH")) %>% filter(qval < FDR)

### 周りのmatrixによるスコアと、control >= 1でもfiltering
D_table <- D_table %>% filter(okay == 1 & control >= T_control) %>% select(-okay)
cat(paste("Filtered number: ", D_table %>% nrow(), "\n"))


#=============================================
# 周りの領域をscanしてピークとなる領域を検出
#=============================================
### 最大のスコアを持つものに入れ替える
D_table2 <- NULL
for(k in 1:nrow(D_table)){
  cen_x <- D_table[k,"id1"]
  cen_y <- D_table[k,"id2"]
  list_x <- max(c(1, cen_x - 5)):min(c(cen_x + 5, nrow(map)))
  list_y <- max(c(1, cen_y - 5)):min(c(cen_y + 5, nrow(map)))
  comb <- expand.grid(list_x, list_y)
  comb <- data.frame(id1=as.integer(comb[,1]), id2=as.integer(comb[,2]), stringsAsFactors = FALSE)
  comb <- dplyr::left_join(comb, D_table, by=c("id1", "id2"))
  comb <- comb %>% filter(!is.na(qval))
  
  ### 少なくとも周りに１つ以上他のsignificantなピークがあるものだけ選ぶ
  if(nrow(comb) > 1){
    comb <- comb %>% filter(qval == min(comb$qval))
    D_table2 <- rbind(D_table2, comb)
  }
}

D_table <- D_table2
rm(D_table2)

### 重複しているentryを除く
D_table <- D_table %>% distinct(id1, id2, .keep_all = TRUE)
cat(paste("After scan surrounded area: ", D_table %>% nrow(), "\n"))
rm(map)


D_table <- dplyr::left_join(D_table, Location %>% rename(chr1=chr, start1=start, end1=end), by=c("id1"="id"))
D_table <- dplyr::left_join(D_table, Location %>% rename(chr2=chr, start2=start, end2=end), by=c("id2"="id"))
D_table <- D_table[,c("chr1", "start1", "end1", "chr2", "start2", "end2", "distance", "score", "control", "fc", "pval", "qval")]
D_table <- D_table %>% arrange(qval)
  
write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






