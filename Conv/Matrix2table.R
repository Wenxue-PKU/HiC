# matrixデータをtableとして出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix files separated by commna(,)"),
  make_option(c("-o", "--out"),help="output file name"),
  make_option(c("--distance"), default="NA", help="max distance, inter for interchromosome"),
  make_option(c("--name"),help="name of sample")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILEs <- unlist(strsplit(as.character(opt["in"]), ","))
NAME_LIST <- unlist(strsplit(as.character(opt["name"]), ","))
FILE_out <- as.character(opt["out"])
FILE_header <- tools::file_path_sans_ext(FILE_out)

if(length(FILEs) != length(NAME_LIST)){
  cat("Number of maps and NAME_LIST are different\n")
  q()
}

SAMPLE_NUMBER <- length(FILEs)

# r.commonの抽出
for(i in 1:SAMPLE_NUMBER){
  FILE_object <- sub(".matrix", ".rds", FILEs[i])
  if(!file.exists(FILE_object)){
    D <- as.matrix(read.table(FILEs[i], header=TRUE, check.names = FALSE))
  }else{
    D <- readRDS(FILE_object)
  }
  r <- rownames(D)
  if(i == 1){
    r.common <- r
  }else{
    r.common <- intersect(r.common, r)
  }
}


# 距離を制限したindexを作成
if(as.character(opt["distance"]) == "inter"){
  index_list <- c()
  for(i in 1:(length(r.common))){
    index <- cbind(r.common[i], r.common[c(i:(length(r.common)))])
    index_list <- rbind(index_list, index)
  }
}else{
  LocList <- strsplit(r.common, ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
  CHROMOSOMES <- unique(as.character(LocMatrix[,1]))
  
  MIN_INDEX = tapply(1:length(r.common), as.character(LocMatrix[,1]), min)
  MAX_INDEX = tapply(1:length(r.common), as.character(LocMatrix[,1]), max)
  THR_LENGTH = tapply(as.numeric(LocMatrix[,3]), as.character(LocMatrix[,1]), length) - 1
  if(as.character(opt["distance"]) != "NA"){
    THR_LENGTH[CHROMOSOMES] = min(as.numeric(as.character(opt["distance"])) / Resolution, THR_LENGTH[CHROMOSOMES])
  }
  index_list <- c()
  for(chr in CHROMOSOMES){
    getindex <- function(d){
      index1 <- MIN_INDEX[chr]:(MAX_INDEX[chr] - d)
      index2 <- index1 + d
      index3 <- cbind(r.common[index1], r.common[index2])
      index3
    }
    tmp <- do.call(rbind, lapply(seq(0, THR_LENGTH[chr]), getindex))
    index_list <- rbind(index_list, tmp)
  }
  rm(tmp)
}
write.table(index_list, file=paste(FILE_header, "_tmp_column.txt", sep=""), sep="\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE, eol = "\n")


# 一旦メモリーを空ける
if(as.character(opt["distance"]) != "inter"){
  rm(LocList, LocMatrix, r.common)
}

# 各ファイルごとにScoreを出力
for(i in 1:SAMPLE_NUMBER){
  FILE_object <- sub(".matrix", ".rds", FILEs[i])
  if(!file.exists(FILE_object)){
    D <- as.matrix(read.table(FILEs[i], header=TRUE, check.names = FALSE))
  }else{
    D <- readRDS(FILE_object)
  }
  
  Score <- as.numeric(D[index_list])
  write.table(Score, file=paste(FILE_header, "_tmp_", NAME_LIST[i], ".txt", sep=""), sep="\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE, eol = "\n")
}


# 全ファイルを結合する
HEADER_list <- c("Loc1", "Loc2", NAME_LIST)
cat(paste(HEADER_list, collapse = "\t"), "\n", file = FILE_out)
FILE_list <- c(paste(FILE_header, "_tmp_column.txt", sep=""), paste(FILE_header, "_tmp_", NAME_LIST, ".txt", sep=""))
system(paste("paste", paste(FILE_list, collapse = " "), ">>", FILE_out))

for(f in FILE_list){
  file.remove(f)
}













