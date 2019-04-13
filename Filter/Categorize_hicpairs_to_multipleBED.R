#!/usr/bin/Rscript
# HiCのペアの情報と、複数のBEDファイルの情報を関連付けるプログラム
# 例:Significantなピークをカテゴリー分けして、さらに、promoter, enhancerなども、overlapしているpeakとの情報を出力


suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-c", "--hic"), default="NA", help="hic pair information"),
  make_option(c("-b", "--BED"), default="NA", help="BED files. Separated by ,"),
  make_option(c("-s", "--sign"), default="NA", help="Characters assign for each BED file. separated by ,. low to high priority order. ex.E,S,P for E < S < P"),
  make_option(c("-o", "--out"), default="NA", help="output directory. add / at the end")
)
opt <- parse_args(OptionParser(option_list=option_list))

# 
# #=============================================
# # test data
# #=============================================
# FILE_hic <- "X:/hideki_projects/2019-01-18_HiC_human_senescent2/out/2019-03-26_significant_association_search_surrounded/IMR90_G_mixAll.txt"
# D_hic <- read.table(FILE_hic, header=TRUE, stringsAsFactors = FALSE, sep="\t", quote="", nrows = 1000)
# FILE_BEDs <- paste0("X:/sharedata/240/HiC_chipseq/", c("Enhancer_G_minus_super.bed", "Super_enhancer_G.bed", "TSS_active_G.bed"))
# SIGNs <- c("E", "S", "P")
# #=============================================



FILE_hic <- as.character(opt["hic"])
FILE_BEDs <- unlist(strsplit(as.character(opt["BED"]), ","))
SIGNs <- unlist(strsplit(as.character(opt["sign"]),","))
DIR_out <- as.character(opt["out"])

if(length(SIGNs) != length(FILE_BEDs)){
  cat("sign and BED file should be same number\n")
  q()
}


suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(dplyr)))

D_hic <- read.table(FILE_hic, header=TRUE, sep="\t", stringsAsFactors = FALSE)
G_L <- GRanges(D_hic %>% dplyr::rename(chr=chr1, start=start1, end=end1) %>% select(chr, start, end))
G_R <- GRanges(D_hic %>% dplyr::rename(chr=chr2, start=start2, end=end2) %>% select(chr, start, end))


con_bin_BED <- function(file, name){
  DATA_head <- read.table(file, header=FALSE, nrows = 5, stringsAsFactors = FALSE, sep="\t", quote="")
  classes <- sapply(DATA_head, class)
  classes[c(4:length(classes))] <- "NULL"
  df <- read.table(file, header=FALSE, colClasses = classes, stringsAsFactors = FALSE, sep="\t", quote="")
  colnames(df) <- c("chr", "start", "end")
  gg <- GRanges(df)
  f1 <- countOverlaps(G_L, gg)
  f2 <- countOverlaps(G_R, gg)
  df <- rbind(data.frame(pairID=which(f1 > 0), side="left", type=name, stringsAsFactors = FALSE),
              data.frame(pairID=which(f2 > 0), side="right", type=name, stringsAsFactors = FALSE))
  df
}

D_signs <- NULL
for(i in 1:length(SIGNs)){
  df <- con_bin_BED(FILE_BEDs[i], SIGNs[i])
  D_signs <- rbind(D_signs, df)
}
rm(df)

D_signs <- D_signs %>% mutate(type=factor(type, levels = c("N", SIGNs), ordered = TRUE))
D_signs.left <- D_signs %>% filter(side=="left") %>% group_by(pairID) %>% summarise(left=max(type)) %>% as.data.frame()
D_signs.right <- D_signs %>% filter(side=="right") %>% group_by(pairID) %>% summarise(right=max(type)) %>% as.data.frame()
rm(D_signs)

D_hic <- dplyr::left_join(D_hic %>% mutate(pairID=row_number()), D_signs.left, by="pairID")
D_hic <- dplyr::left_join(D_hic, D_signs.right, by="pairID")
rm(D_signs.left, D_signs.right)
D_hic <- D_hic %>% mutate(left=ifelse(is.na(left), "N", as.character(left)), right=ifelse(is.na(right), "N", as.character(right)))
D_hic <- D_hic %>% mutate(left=factor(left, levels = c("N", SIGNs), ordered = TRUE), right=factor(right, levels = c("N", SIGNs), ordered = TRUE))
D_hic <- D_hic %>% mutate(pattern=ifelse(left > right, paste(as.character(right), as.character(left), sep="-"), 
                                         paste(as.character(left), as.character(right), sep ="-")))

write.table(D_hic %>% select(-pairID), paste0(DIR_out, "pairs.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


pattern_comb <- expand.grid(c("N", SIGNs), c("N", SIGNs))
pattern_comb <- rev(paste(pattern_comb[,1], pattern_comb[,2], sep="-"))

D_out <- D_hic %>% mutate(pattern=factor(pattern, levels=pattern_comb)) %>% arrange(pattern)%>% group_by(pattern) %>% tally() %>% as.data.frame()
write.table(D_out, paste0(DIR_out, "Category_peaks.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



### 各BEDファイルについて、対応するカテゴリーを出力する
output_category_for_BED <- function(file, name){
  D_bed <- read.table(file, header=FALSE, stringsAsFactors = FALSE, sep="\t", quote="")
  colnames(D_bed)[1:3] <- c("chr", "start", "end")
  gg <- GRanges(D_bed[,1:3])
  f1 <- findOverlaps(gg, G_L)
  f2 <- findOverlaps(gg, G_R)
  D_bed <- D_bed %>% mutate(bedID=row_number())
  df <- rbind(data.frame(bedID=f1@from, pattern=D_hic$pattern[f1@to], stringsAsFactors = FALSE),
              data.frame(bedID=f2@from, pattern=D_hic$pattern[f2@to], stringsAsFactors = FALSE))
  df <- df %>% mutate(pattern=factor(pattern, levels=pattern_comb, ordered = TRUE))
  
  ### 数のまとめ
  D_out.allowOverlap <- df %>% group_by(pattern) %>% summarize(allowOverlap=n()) %>% as.data.frame()
  D_out.unique <- df %>% group_by(bedID) %>% summarize(pattern=max(pattern)) %>% as.data.frame() %>% group_by(pattern) %>% summarize(unique=n()) %>% as.data.frame()
  D_out <- dplyr::left_join(D_out.allowOverlap, D_out.unique, by="pattern")
  D_out <- D_out %>% mutate(pattern=as.character(pattern))
  D_out <- rbind(D_out, c("Total overlapped", colSums(D_out[,-1])))
  
  df <- df %>% group_by(bedID) %>% summarize(pattern=paste(unique(pattern), collapse = ",")) %>% as.data.frame()
  D_bed <- dplyr::left_join(D_bed, df, by="bedID")
  write.table(D_bed %>% select(-bedID), paste0(DIR_out, "BED_", name, ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "")
  
  D_out <- rbind(D_out, c("Total not overlapped", rep(sum(is.na(D_bed$pattern)) ,2)))
  write.table(D_out, paste0(DIR_out, "Category_", name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

for(i in 1:length(SIGNs)){
  output_category_for_BED(FILE_BEDs[i], SIGNs[i])
}




