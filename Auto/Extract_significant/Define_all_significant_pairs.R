#!/usr/bin/Rscript
# Extract significant associations


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output file"),
  make_option(c("--min"), default="20000", help="minimum distance to check(default : 20kb)"),
  make_option(c("--max"), default="2000000", help="maximum distance to check (default : 2Mb)")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))

FILE_out <- as.character(opt["out"])
FILE_matrix <- as.character(opt["in"])
T_max <- as.numeric(as.character(opt["max"]))
T_min <- as.numeric(as.character(opt["min"]))
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


SEARCH_AREA <- 5


#=============================================
# 周りのスコアを見てフィルタリング
#=============================================
mask <- data.frame(
          V1 = c(1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 2),
          V2 = c(1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 2),
          V3 = c(1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 2),
          V4 = c(1, 1, 1, 0, 0, 0, 0, 0, 2, 2, 2),
          V5 = c(1, 1, 1, 0, -1, -1, -1, 0, 2, 2, 2),
          V6 = c(0, 0, 0, 0, -1, -2, -1, 0, 0, 0, 0),
          V7 = c(1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1),
          V8 = c(1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1),
          V9 = c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1),
         V10 = c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1),
         V11 = c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1)
)
mask <- as.matrix(mask)
mask_all <- which(mask>-2, arr.ind = TRUE)
mask_center <- which(mask<0, arr.ind = TRUE)
mask_surround <- which(mask>0, arr.ind = TRUE)
mask_leftBottom <- which(mask==2, arr.ind = TRUE)
checkSurrounded <- function(i,j){
  ### -6するのは中心をずらすため
  index_all <- cbind(mask_all[,1]+i-6, mask_all[,2]+j-6)
  index_center <- cbind(mask_center[,1]+i-6, mask_center[,2]+j-6)
  index_surround <- cbind(mask_surround[,1]+i-6, mask_surround[,2]+j-6)
  index_leftBottom <- cbind(mask_leftBottom[,1]+i-6, mask_leftBottom[,2]+j-6)
  
  sss_middle <- as.numeric(map[i,j])
  sss_center <- as.numeric(map[index_center]) 
  
  if(is.na(sss_middle)){
    0
  }else{
    ### 90% 以上がスコアで埋まっている
    sss <- as.numeric(map[index_all])
    FLAG_fill <- (sum(!is.na(sss) & sss > 0) / length(sss)) > 0.90
    
    ### centerは少なくとも4つ以上
    FLGA_center <- sum(!is.na(sss_center)) > 3

    if(!(FLAG_fill && FLGA_center)){
      0
    }else{
      ### 周りのスコアよりも中央付近のスコアが高い
      sss_surround <- as.numeric(map[index_surround])
      FLAG_surround <- mean(sss_surround, na.rm = TRUE) > mean(sss_center, na.rm = TRUE)
      if(FLAG_surround){
        0
      }else{
        ### Left bottom cornerよりも高い
        sss_leftBottom <- as.numeric(map[index_leftBottom])
        FLAG_leftBottom <- mean(sss_leftBottom, na.rm = TRUE) < sss_middle
        if(FLAG_leftBottom){
          1
        }else{
          0
        }
      }
    }
  }
}


D_table <- c()
MAX_NUM <- min(c(NUM_LINE-1, as.integer(T_max/Resolution)))
MIN_NUM <- max(c(1, as.integer(T_min/Resolution)))
for(d in MIN_NUM:MAX_NUM){
  index1 <- (SEARCH_AREA + 1):(NUM_LINE - d - SEARCH_AREA)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  
  wrapper <- function(k){
    checkSurrounded(as.integer(index3[k,1]), as.integer(index3[k,2]))
  }
  Okay <- sapply(1:nrow(index3), wrapper)
  
  scores <- as.numeric(map[index3])
  Average <- mean(scores, na.rm=TRUE)
  if(!is.na(Average) & Average != 0){
    distance <- d*Resolution
    Pvalues <- pnorm(scores, mean=Average, sd=sd(scores, na.rm = TRUE), lower.tail = FALSE)
    df <- data.frame(id1=index1, id2=index2, distance=distance, score=scores, control=Average, 
                     fc=log2(scores/Average), pval=Pvalues, okay=Okay, stringsAsFactors = FALSE)
    df <- df %>% filter(!is.na(Pvalues))
    D_table <- rbind(D_table, df)
  }
}
rm(df, map, index3)

NUM_all <- D_table %>% nrow()
D_table <- D_table %>% mutate(qval=p.adjust(pval, method = "BH")) %>% filter(qval < 0.01)

### 周りのmatrixによるスコアと、control >= 1というでもfiltering
D_table <- D_table %>% filter(okay == 1 & control >= 1) %>% select(-okay)

NUM_sig <- D_table %>% nrow()
cat(paste(NUM_sig, NUM_all, sep="\t"))

D_table <- dplyr::left_join(D_table, Location %>% rename(chr1=chr, start1=start, end1=end), by=c("id1"="id"))
D_table <- dplyr::left_join(D_table, Location %>% rename(chr2=chr, start2=start, end2=end), by=c("id2"="id"))
D_table <- D_table[,c("chr1", "start1", "end1", "chr2", "start2", "end2", "distance", "score", "control", "fc", "pval", "qval")]
D_table <- D_table %>% arrange(qval)
  
write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)






