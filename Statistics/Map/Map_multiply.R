# Map全体に一定の値をかけたものを出力する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), help="input matrix file"),
  make_option(c("-o", "--out"), help="output matrix file or rds file "),
  make_option(c("-v", "--value"), default="1", help="value to multiply")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_matrix <- as.character(opt["in"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
if(!file.exists(FILE_object)){
  map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
}else{
  map <- readRDS(FILE_object)
}

VALUE_multiply <- as.numeric(as.character(opt["value"]))
map <- map * VALUE_multiply


FILE_out <- as.character(opt["out"])

if(sum((grep("\\.rds$", FILE_out))) == 1){
  saveRDS(map, FILE_out)
}else{
  cat("\t", file=FILE_out)
  suppressWarnings( write.table(map, file=FILE_out, append=TRUE, quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=TRUE))
}
