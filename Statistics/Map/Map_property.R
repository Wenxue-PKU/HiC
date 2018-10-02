# Mapの統計データを出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), help="input matrix file separated by ,"),
  make_option(c("-o", "--out"), default="NA", help="output histogram"),
  make_option(c("--lineMAX"), default="FALSE", help="max line score instead of individual"),
  make_option(c("--lineSUM"), default="FALSE", help="sum line score instead of individual"),
  make_option(c("--normalize"), default="FALSE", help="distance normalize or not"),
  make_option(c("--distance"), default="NA", help="specific distance"),
  make_option(c("--info"), default="NA", help="what kind of information will be output")
)
opt <- parse_args(OptionParser(option_list=option_list))



DistanceNormalize <- function(map){
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
  map
}

Filter_distance <- function(map){
  r <- rownames(map)
  LocList <- strsplit(r, ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  rownames(LocMatrix) <- r
  
  comb <- expand.grid(r, r)
  Distance <- abs(as.numeric(LocMatrix[comb[,1],2]) - as.numeric(LocMatrix[comb[,2],2]))
  index_NG <- Distance != as.numeric(as.character(opt["distance"]))
  map[cbind(comb[index_NG,1], comb[index_NG,2])] <- NA
  map
}


FILEs <- strsplit(as.character(opt["in"]), ",")

Score <- c()
for(FILE_matrix in unlist(strsplit(as.character(opt["in"]), ","))){
  FILE_object <- sub(".matrix", ".rds", as.character(FILE_matrix))
  if(!file.exists(FILE_object)){
    map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  }else{
    map <- readRDS(FILE_object)
  }
  if(as.character(opt["normalize"]) == "TRUE"){
    map <- DistanceNormalize(map)
  }
  if(as.character(opt["distance"]) != "NA"){
    map <- Filter_distance(map)
  }
  if(as.character(opt["lineMAX"]) == "TRUE"){
    Score <- c(Score, apply(map, 1, max, na.rm=TRUE))
  }else if(as.character(opt["lineSUM"]) == "TRUE"){
    Score <- c(Score, apply(map, 1, sum, na.rm=TRUE))
  }else{
    Score <- c(Score, as.numeric(map[upper.tri(map, diag=TRUE)]))
  }
}

Score <- Score[is.finite(Score)]
Score <- Score[!is.na(Score)]
Score <- round(Score)

FILE_out <- as.character(opt["out"])
if(FILE_out != "NA"){
  b <- seq(-1, max(Score)+1, by=1)
  t <- hist(round(Score), breaks=b, plot=FALSE, include.lowest = FALSE)
  OUT <- cbind(t$breaks[-1], t$counts)
  colnames(OUT) <- c("Read", "Count")
  write.table(OUT, file=FILE_out, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE, eol = "\n")
}


if(as.character(opt["info"]) != "NA"){
  switch(as.character(opt["info"]),
         Summary={
           cat("SUM:\t", sum(Score), "\n", sep="")
           cat("Average:\t", mean(Score), "\n", sep="")
           cat("Max:\t", max(Score), "\n", sep="")
           cat("Min:\t", min(Score), "\n", sep="")
           cat("% of 0:\t", sum(Score==0)/length(Score)*100, "\n", sep="")
         },
         SUM={
           cat(sum(Score),"\n")
         },
         Average={
           cat(mean(Score), "\n")
         },
         MAX={
           cat(max(Score), "\n")
         },
         MIN={
           cat(min(Score), "\n")
         },
         SD={
           cat(sd(Score), "\n")
         },
         SDAve={
           cat(sd(Score) / mean(Score), "\n")
         },
         zero={
           cat(sum(Score==0)/length(Score)*100, "\n")
         },
         {
           cat("Please select info from Summary, SUM, Average, MAX, MIN, SD, SDAve, zero\n")
         }
  )
}






