# 顕著に変化している組合せを抽出する(100kb - 3Mbの範囲)

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--file1"),help="matrix file1"),
  make_option(c("-b", "--file2"),help="matrix file2"),
  make_option(c("-p", "--pval"), defaul=0.01, help="p-value threshold"),
  make_option(c("-f", "--fold"), defaul=1.5, help="fold-change threshold"),
  make_option(c("--min"), defaul=100000, help="minimum distance to check"),
  make_option(c("--max"), defaul=3000000, help="maximum distance to check"),
  make_option(c("-v", "--vervose"), default="FALSE", help="output intermediate status"),
  make_option(c("-o", "--out"),help="output text file")
)
opt <- parse_args(OptionParser(option_list=option_list))

FLAG_status <- eval(parse(text=as.character(opt["vervose"])))
FILE_output <- as.character(opt["out"])

FILE_matrix1 <- as.character(opt["file1"])
FILE_matrix2 <- as.character(opt["file2"])
FILE_object1 <- sub("matrix", "rds", FILE_matrix1)
FILE_object2 <- sub("matrix", "rds", FILE_matrix2)

# load map1
if(FLAG_status){
  cat("load map1 ... ")
}
if(!file.exists(FILE_object1)){
  map1 <- as.matrix(read.table(FILE_matrix1, header=TRUE, check.names = FALSE))
}else{
  map1 <- readRDS(FILE_object1)
}
r <- rownames(map1)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[2,2]) - as.numeric(LocMatrix[1,2])
if(FLAG_status){
  cat(Resolution/1000, "kb resolution, Length", max(as.numeric(LocMatrix[,2])), "\n")
}
rm(LocList,LocMatrix)



# 組合せを作成する
if(FLAG_status){
  cat("make combinations of ", as.character(opt["min"]), " < distance < ", as.character(opt["max"]), " ... ")
}
index_target <- c()
distance <- c()
for(d in 1:(nrow(map1)-1)){
  ds<- d * Resolution
  if(ds > as.numeric(as.character(opt["max"]))){
    break
  }
  if(ds> as.numeric(as.character(opt["min"]))){
    index_1 <- 1:(nrow(map1) - d)
    index_2 <- index_1 + d
    index_3 <- cbind(r[index_1], r[index_2])
    index_target <- rbind(index_target, index_3)
    distance <- c(distance, rep(ds, length(index_1)))
  }
}
if(FLAG_status){
  cat(length(distance), " combinations\n")
}


# load map2
if(FLAG_status){
  cat("load map2 ... ")
}
if(!file.exists(FILE_object2)){
  map2 <- as.matrix(read.table(FILE_matrix2, header=TRUE, check.names = FALSE))
}else{
  map2 <- readRDS(FILE_object2)
}
r2 <- rownames(map2)
ok <- index_target[,1] %in% r2 & index_target[,2] %in% r2
index_target <- index_target[ok,]
distance <- distance[ok]
rm(r2, ok)
map1 <- map1[index_target]
map2 <- map2[index_target]
if(FLAG_status){
  cat(length(map2), " targets\n")
}


# ratioの計算
if(FLAG_status){
  cat("Calc ratio ... ")
}
ratio <- log2(map1 / map2)
no <- is.na(ratio) | is.infinite(ratio)
index_target <- index_target[!no,]
ratio <- ratio[!no]
distance <- distance[!no]
map1 <- map1[!no]
map2 <- map2[!no]
if(FLAG_status){
  cat(length(distance), " combinations\n")
}
fd <- 2**ratio
fd <- ifelse(ratio < 0, 1/fd, fd)

# 予めbackgroundの平均値と標準偏差を計算しておく
if(FLAG_status){
  cat("prepare background ... ")
}
bg_dis <- sort(unique(distance))
getMean <- function(i){
  mean(ratio[distance == bg_dis[i]])
}
getSd <- function(i){
  sd(ratio[distance == bg_dis[i]])
}
bg_mean <- sapply(1:length(bg_dis), getMean)
bg_sd <- sapply(1:length(bg_dis), getSd)
if(FLAG_status){
  cat("ok\n")
}




# P-valueを計算する
if(FLAG_status){
  cat("pvalue calculation ... ")
}
getPval <- function(i){
  n <- which(bg_dis == distance[i])
  p <- pnorm(ratio[i], mean=bg_mean[n], sd=bg_sd[n])
  if(is.na(p)){
    p <- 1
  }
  p
  # cat(ratio[i], p, bg_mean[n], bg_sd[n], "\n", sep = "\t")
}
pval <- sapply(1:length(ratio), getPval)
pval <- ifelse(pval > 0.5, 1-pval, pval)
if(FLAG_status){
  cat("ok\n")
}


# significantを指定されたp-value以下でmap1, map2のスコアが0以上、ratioが指定された値以上とする
T_pval <- as.numeric(opt["pval"])
T_fd <- as.numeric(opt["fold"])
if(FLAG_status){
  cat("Significant were defind as p-value < ", T_pval, " and fold > ", T_fd, " and map1>0 map2>0\n", sep="")
}

sig <- pval < T_pval & fd > T_fd & map1 > 0 & map2 > 0
if(FLAG_status){
  cat(sum(sig), "/", length(sig), " is significant (", round(sum(sig)/length(sig)*100, digits = 3), "%)\n")
}
if(sum(sig) == 0){
  q()
}

# 出力
if(FLAG_status){
  cat("output result as", FILE_output, " ... ")
}


OUT <- cbind(index_target[sig,1], index_target[sig,2], paste(distance[sig]/1000, "kb", sep=""), format(map1[sig], digits=3), format(map2[sig], digits=3),
             round(ratio[sig], digits = 3), round(fd[sig], digits = 3), format(pval[sig], digits = 3))
colnames(OUT) <- c("loc1", "loc2", "distance", "map1", "map2",  "log2(map1/map2)", "map1/map2", "pval")
OUT <- OUT[order(pval[sig], 1-abs(ratio[sig])),]
write.table(OUT, FILE_output, col.names = TRUE, row.names = FALSE, append = FALSE, sep = "\t", quote = FALSE, eol = "\n")
if(FLAG_status){
  cat("finished\n")
}






