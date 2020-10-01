#!/usr/bin/Rscript
# PCA scoreの計算

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file (rds or text)"),
  make_option(c("--format"), default="rds", help="input file format (default: rds)"),
  make_option(c("-o", "--out"), default="NA", help="output file"),
  make_option(c("-g", "--gene"), default="NA", help="gene bed file")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(GenomicRanges)))


options(scipen=10)
FILE_out <- as.character(opt["out"])
FILE_matrix <- as.character(opt["in"])
FILE_format <- as.character(opt["format"])


FILE_gene <- as.character(opt["gene"])
if(FILE_gene== "NA"){
  cat("gene bed file is reuiqred")
  q()
}

if(FILE_format == "rds"){
  map <- readRDS(FILE_matrix)
}else if(FILE_format == "matrix"){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  cat("Unknown input file format")
  q()
}
r <- rownames(map)
D_col <- as.data.frame(r) %>% mutate(name=r) %>% tidyr::separate(r, c("chr", "start", "end"), ":")
D_col <- D_col %>% mutate(start=as.numeric(start), end=as.numeric(end))
G_col <- GRanges(D_col)


# 遺伝子のcount
D_gene <- fread(FILE_gene, header=F)
D_gene <- D_gene[,1:3]
colnames(D_gene) <- c("chr", "start", "end")
G_gene <- GRanges(D_gene)
ov <- countOverlaps(G_col, G_gene)
D_col <- D_col %>% mutate(geneNumber=ov)
rm(G_col, G_gene, D_gene)

RESOLUTION <- D_col %>% mutate(resolution=end-start+1) %>% head(n=1) %>% pull(resolution)


# Observed / Expectのmatrixに変換する
map_expect <- map
NUM_LINE <- nrow(map)
for(d in 0:(NUM_LINE-1)){
  index1 <- 1:(NUM_LINE - d)
  index2 <- index1 + d
  index3 <- cbind(index1, index2)
  index4 <- cbind(index2, index1)
  Average <- mean(as.numeric(map[index3]), na.rm=TRUE)
  map_expect[index3] <- Average
  map_expect[index4] <- Average
}
map <- ifelse(map_expect == 0, 1, map / map_expect)

# もしsdが0だったらNAにする
sdlist <- apply(map,1,sd, na.rm=TRUE)
index <- which(sdlist == 0)
map[index, ] <- NA
map[, index] <- NA

getCor <- function(i){
  apply(map, 2, function(x) { cor(x, map[,i], use="pairwise.complete.obs", method="pearson")})
}
options(warn=-1)
map <- rbind(sapply(1:ncol(map), getCor))
colnames(map) <- rownames(map)
options(warn=1)

SUM <- apply(is.na(map), 2, sum)
NonNABin <- r[SUM != nrow(map)]

# PCA analysis
pca <- prcomp(map[NonNABin, NonNABin],  scale=TRUE)


# 合計で何個の連続したcompartmentがあるのかを数え、もし、その数が10未満で、第2主成分では、10以上であった場合、第2主成分を用いる
pca_score <- rep(NA, nrow(map))
FirstComponent <- pca$x[,1]
SecondComponent <- pca$x[,2]
names(pca_score) <- r

checkCompartmentNum <- function(scores){
  scores <- sign(scores)
  scores <- ifelse(is.na(scores), 0, scores)
  NUM_total_compartment <- 1
  old_sign <- scores[1]
  for(i in 2:length(scores)){
    new_sign <- scores[i]
    if(new_sign != old_sign){
      NUM_total_compartment  <- NUM_total_compartment + 1
    }
    old_sign <- new_sign
  }
  NUM_total_compartment
}
if(checkCompartmentNum(FirstComponent) < 10 && checkCompartmentNum(SecondComponent) > 10){
  pca_score[NonNABin] <- SecondComponent
}else{
  pca_score[NonNABin] <- FirstComponent
}


# 遺伝子の密度が多い方をCompartmentA、低い方をcompartmentBとする
group1 <- intersect(names(which(pca_score > 0)), r)
group2 <- intersect(names(which(pca_score < 0)), r)
GeneAve_group1 <- D_col %>% filter(name %in% group1) %>% pull(geneNumber) %>% mean()
GeneAve_group2 <- D_col %>% filter(name %in% group2) %>% pull(geneNumber) %>% mean()

if(GeneAve_group1 < GeneAve_group2){
  # group2がcompartmentAなのでスコアを正に調整
  if(sum(pca_score[group2], na.rm=TRUE) < 0){
    pca_score <- pca_score * (-1)
  }
}else{
  # group1がcompartmentAなのでスコアを正に調整
  if(sum(pca_score[group1], na.rm=TRUE) < 0){
    pca_score <- pca_score * (-1)
  }
}

D_out <- data.frame(name=names(pca_score), pca=pca_score) %>% tidyr::separate(name, c("chr", "start", "end"), ":")
D_out <- D_out %>% mutate(comp=if_else(pca>0, "A", "B"))
write.table(D_out, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



