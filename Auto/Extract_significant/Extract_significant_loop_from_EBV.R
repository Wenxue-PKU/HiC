#!/usr/bin/Rscript
# significantなピークを定義する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="input matrix"),
  make_option(c("-o", "--out"), default="NA", help="significant peak list")
)
opt <- parse_args(OptionParser(option_list=option_list))



suppressWarnings(suppressMessages(library(dplyr)))

FILE_matrix <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])

map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
r <- rownames(map)
r2id <- 1:length(r)
names(r2id) <- r
NUM_LINE <- nrow(map)

splitLocation <- function(r){
  LocList <- strsplit(r, ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  Location <- data.frame(chr=LocMatrix[,1], start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]), stringsAsFactors = FALSE)
  Location
}
D_mat <- as.data.frame(map) %>% mutate(bin1=r) %>% tidyr::gather(key=bin2, value=score, -bin1)
Loc_L <- splitLocation(D_mat$bin1)
Loc_R <- splitLocation(D_mat$bin2)
D_mat <- D_mat %>% mutate(id1=r2id[D_mat$bin1], start1=Loc_L$start, end1=Loc_L$end, id2=r2id[D_mat$bin2], start2=Loc_R$start, end2=Loc_R$end)
Sliding_window <- Loc_L[2,"end"] - Loc_L[1,"end"]
rm(Loc_L, Loc_R)
D_loc <- D_mat %>% distinct(id1, id2, .keep_all = TRUE) %>% select(id1, id2, start1, end1, start2, end2)
D_mat <- D_mat %>% mutate(distance=abs(start1-start2))




### NAを除く
D_mat <- D_mat %>% filter(!is.na(score))
D_mat <- D_mat %>% filter(start1 < start2)
D_mat <- D_mat %>% mutate(dcate=dplyr::ntile(distance, 40)) 



MAX_NUM <- NUM_LINE - 1
MIN_NUM <- 1
D_table <- c()
for(d in MIN_NUM:MAX_NUM){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  
  scores <- as.numeric(map[index3])
  
  ### background score
  ### top1%のデータだけbackgroundから取り除く
  distance_for_dmat <- D_mat %>% filter(distance==d*Sliding_window) %>% pull(dcate) %>% unique()
  scores_back <- D_mat %>% filter(dcate %in% distance_for_dmat) %>% pull(score)
  scores_back <- scores_back[scores_back < quantile(scores_back, prob=0.99)]
  distance_category <- paste(D_mat %>% filter(dcate %in% distance_for_dmat) %>% pull(distance) %>% min(na.rm=TRUE),
                             D_mat %>% filter(dcate %in% distance_for_dmat) %>% pull(distance) %>% max(na.rm=TRUE),
                             sep="-")
  
  Average <- mean(scores_back, na.rm = TRUE)
  if(!is.na(Average) & Average != 0){
    distance <- d*Sliding_window
    Pvalues <- pnorm(scores, mean=Average, sd=sd(scores_back), lower.tail = FALSE)
    df <- data.frame(id1=index1, id2=index2, distance=distance, distance_category, score=scores, control=Average, 
                     fc=log2(scores/Average), pval=Pvalues, stringsAsFactors = FALSE)
    df <- df %>% filter(!is.na(Pvalues))
    D_table <- rbind(D_table, df)
  }
}
rm(df, index1, index2, index3, scores, Pvalues)

D_table <- D_table %>% mutate(qval=p.adjust(pval, method = "BH"))
D_table <- D_table %>% filter(fc > 0.5)

D_out <- dplyr::left_join(D_table, D_loc, by=c("id1", "id2"))


D_out <- D_out[,c("start1", "end1",  "start2", "end2", "distance", "distance_category", "score", "control", "fc", "pval", "qval")]
D_out <- D_out %>% arrange(qval)

write.table(D_out, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

