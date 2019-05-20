#!/usr/bin/Rscript
# MatrixのMoving averageを計算

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), help="original matrix file"),
  make_option(c("-o", "--out"), help="output matrix file or rds file "),
  make_option(c("--keep_na"), default="FALSE", help="convert NA to 0 (FALSE) or keep (TRUE)"),
  make_option(c("-n", "--distance"), default=1, help="number of surrounded bins to calculate average")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_matrix <- as.character(opt["in"])
FILE_object <- sub("matrix", "rds", FILE_matrix)
FILE_out <- as.character(opt["out"])
dd <- as.numeric(as.character(opt["distance"]))


if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}

if(!eval(parse(text=as.character(opt["keep_na"])))){
  map[is.na(map)] <- 0
}

map_new <- map
for (i in (1 : (nrow(map)))){
  for (j in (1 : (ncol(map)))){
    map_new[i,j] = mean(map[max(1,i-dd):min(nrow(map),i+dd),
                               max(1,j-dd):min(ncol(map),j+dd)], na.rm = TRUE)
  }
}

if(sum((grep("\\.rds$", FILE_out))) == 1){
  saveRDS(map_new, FILE_out)
}else{
  suppressWarnings( write.table(map_new, file=FILE_out, quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA))
}


