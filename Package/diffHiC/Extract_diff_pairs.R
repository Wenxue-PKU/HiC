#!/usr/bin/Rscript
# diff HiC

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="all rds separated by ,"),
  make_option(c("-o", "--out"), default="NA", help="output file of p-values"),
  make_option(c("-g", "--group"), default="NA", help="groups (1 or 2 or 3) separated by ,. 3 will be used for estimating dispersion. logFC >0 means 1<2")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(diffHic)))
suppressWarnings(suppressMessages(library(edgeR)))

FILE_samples <- as.character(opt["in"])
FILE_samples <- unlist(strsplit(FILE_samples, ","))
group <- as.character(opt["group"])
group <- unlist(strsplit(group, ","))
FILE_output <- as.character(opt["out"])


if(file.exists(FILE_output)){
  file.remove(FILE_output)
}

SAMPLES <- 1:length(FILE_samples)

getMap <- function(FILE_matrix){
  FILE_object <- sub(".matrix", ".rds", FILE_matrix)
  if(!file.exists(FILE_object)){
    map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  }else{
    map <- readRDS(FILE_object)
  }
  map
}

# cat("Load maps...\n")
r.common <- c()
maps <- list()
for(ss in SAMPLES){
  maps[[ss]] <- getMap(FILE_samples[ss])
  if(length(r.common) == 0){
    r.common <- rownames(maps[[ss]])
  }else{
    r.common <- intersect(r.common, rownames(maps[[ss]]))
  }
}
LocList <- strsplit(r.common, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
df <- data.frame(chr=as.character(LocMatrix[,1]), start=as.numeric(LocMatrix[,2]), end=as.numeric(LocMatrix[,3]))
regions <- makeGRangesFromDataFrame(df)
rm(LocList, df)

# cat("Filter maps...\n")
lib.size <- c()
cm <- list()
for(ss in SAMPLES){
  map.filter <- maps[[ss]][r.common, r.common]
  cm[[ss]] <- ContactMatrix(map.filter, regions, regions)
  lib.size <- c(lib.size, sum(map.filter[upper.tri(map.filter, diag = TRUE)], na.rm = TRUE))
}
rm(maps)


to.keep <- as.matrix(cm[[1]])
for(i in 2:length(cm)){
  to.keep <- to.keep + as.matrix(cm[[i]])
}
to.keep <- to.keep!=0
# cat("using n=", sum(as.numeric(to.keep)>0), " / ", length(as.numeric(to.keep)), "\n")


for(i in 1:length(cm)){
  iset <- deflate(cm[[i]], extract=to.keep)
  if(i==1){
    data <- iset
  }else{
    data <- cbind(data, iset)
  }
}

# cat("Make data ...\n")
interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
data$totals <- lib.size

# cat("edgeR analysis ...\n")
colnames(data) <- SAMPLES
assay(data, "counts") <- assay(data)

keep <- aveLogCPM(asDGEList(data, group=group)) > 0
data <- data[keep,]

# cat("Normalization factor ...\n")
# y <- calcNormFactors(data, type="loess", se.out=TRUE)
y <- calcNormFactors(data)

# design matrixの作成
design <- model.matrix(~factor(group))


# cat("Estimate dispersion...\n")
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)


result <- glmQLFTest(fit, coef=2)
adj.p <- p.adjust(result$table$PValue, method="BH")
NUM <- sum(result$table$PValue <= 0.05)

useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
inter.frame <- as.data.frame(interactions(data))[,useful.cols]
results.r <- data.frame(inter.frame, result$table, FDR=adj.p)
results.r <- results.r[order(results.r$PValue),]
write.table(results.r, file=FILE_output, sep="\t", quote=FALSE, row.names=FALSE, col.names = TRUE, append = FALSE)


