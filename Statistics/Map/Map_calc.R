# matrixの単純な計算

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--file1"), help="matrix file1"),
  make_option(c("-b", "--file2"), help="matrix file2"),
  make_option(c("-o", "--out"), help="output matrix file or rds file "),
  make_option(c("--keep_na"), default="FALSE", help="convert NA to 0 (FALSE) or keep (TRUE)"),
  make_option(c("--method"), default="subtract", help="how to calculate. subtract, add")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_matrix1 <- as.character(opt["file1"])
FILE_matrix2 <- as.character(opt["file2"])
FILE_object1 <- sub("matrix", "rds", FILE_matrix1)
FILE_object2 <- sub("matrix", "rds", FILE_matrix2)


if(!file.exists(FILE_object1)){
  map1 <- as.matrix(read.table(FILE_matrix1, header=TRUE, check.names = FALSE))
}else{
  map1 <- readRDS(FILE_object1)
}
if(!file.exists(FILE_object2)){
  map2 <- as.matrix(read.table(FILE_matrix2, header=TRUE, check.names = FALSE))
}else{
  map2 <- readRDS(FILE_object2)
}

if(!eval(parse(text=as.character(opt["keep_na"])))){
  map1[is.na(map1)] <- 0
  map2[is.na(map2)] <- 0
}

r1 <- rownames(map1)
r2 <- rownames(map2)
r.common <- intersect(r1, r2)
map1 <- map1[r.common, r.common]
map2 <- map2[r.common, r.common]

method <- as.character(opt["method"])
switch(method,
       "subtract" = map <- map1 - map2,
       "add" = map <- map1 + map2,
       cat(method, " is not defined")
)


FILE_out <- as.character(opt["out"])

if(sum((grep("\\.rds$", FILE_out))) == 1){
  saveRDS(map, FILE_out)
}else{
  suppressWarnings( write.table(map, file=FILE_out, quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA))
}


