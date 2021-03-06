#!/usr/bin/Rscript
# Compare two distribution by drawinwg histogram

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="score files, separated by commma"),
  make_option(c("--control"), defaul="NA", help="score files for control"),
  make_option(c("--legend_loc"), default="topright", help="legend location"),
  make_option(c("--column"), default=1, help="column number of score"),
  make_option(c("--colors"), default="NA", help="color of each data"),
  make_option(c("--names"), help="name of each files"),
  make_option(c("--bin_width"), defaul="NA", help="width of histogram. NA split bins to 30. round make round number"),
  make_option(c("--min"), defaul="NA", help="minimum of x-axis"),
  make_option(c("--max"), defaul="NA", help="maxmum of x-axis"),
  make_option(c("--title"), defaul="", help="title of graph"),
  make_option(c("--xlab"), defaul="Read#", help="xlabel"),
  make_option(c("-o", "--out"),help="output file"),
  make_option(c("--sampling"), default="NA", help="NA for not do random sampling. Specify the sampling (repeating) number (ex.10000)"),
  make_option(c("--combination"), default=100, help="if do random sampling. Number of target (ex. 100)")
)
opt <- parse_args(OptionParser(option_list=option_list))


FILE_ins <- c("E-E_SCIMR.txt", "CE-E_SCIMR.txt")
DIR <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/out/2018-06-18_combination_of_10kb_windows2/hic_score/"
FILE_ins <- paste(DIR, FILE_ins, sep="")
NAMEs <- c("target", "control")

FILE_ins <- unlist(strsplit(as.character(opt["in"]), ","))
if(as.character(opt["control"]) != "NA"){
  FILE_controls <- unlist(strsplit(as.character(opt["control"]), ","))
}

if(as.character(opt["colors"]) != "NA"){
  Colors <- unlist(strsplit(as.character(opt["colors"]), ","))
}else{
  Colors <- c('#fc8d62', '#66c2a5','#8da0cb','#e78ac3','#a6d854','#ffd92f')[1:length(FILE_ins)]
}


NAMEs <- unlist(strsplit(as.character(opt["names"]), ","))
FILE_out <- as.character(opt["out"])
COLUMN_NUM <- as.numeric(as.character(opt["column"]))
LEGEND_LOC <- as.character(opt["legend_loc"])
NUM_sampling <- as.character(opt["sampling"])
NUM_combination <- as.numeric(as.character(opt["combination"]))

getSampling <- function(score){
  sapply(1:as.numeric(NUM_sampling), function(i){ mean(sample(score, NUM_combination, replace = FALSE), na.rm=TRUE)})
}


Limit <- list()
DATA <- list()
for(i in 1:length(FILE_ins)){
  D <-  read.table(FILE_ins[i], header=F, stringsAsFactors = FALSE)
  Score <- D[,COLUMN_NUM]
  if(NUM_sampling != "NA"){
    Score <- getSampling(Score)
  }
  if(as.character(opt["control"]) != "NA"){
    D <-  read.table(FILE_controls[i], header=F, stringsAsFactors = FALSE)
    Score_control <- D[,COLUMN_NUM]
    if(NUM_sampling != "NA"){
      Score_control <- getSampling(Score_control)
    }
    Score <- log2((Score + 1)/(Score_control+1))
  }
  DATA[[NAMEs[i]]] <- Score[!is.na(Score)]
  
  if(i==1){
    Limit[['max']] <- max(Score, na.rm = TRUE)
    Limit[['min']] <- min(Score, na.rm = TRUE)
  }else{
    if(max(Score, na.rm = TRUE) > Limit[['max']]){
      Limit[['max']] <- max(Score, na.rm = TRUE)
    }
    if(min(Score, na.rm = TRUE) < Limit[['min']]){
      Limit[['min']] <- min(Score, na.rm = TRUE)
    }
  }
}

if(as.character(opt["max"]) != "NA"){
  MAX <- as.numeric(as.character(opt["max"]))
}else{
  MAX <- round(Limit[['max']])
}
if(as.character(opt["min"]) != "NA"){
  MIN <- as.numeric(as.character(opt["min"]))
}else{
  MIN <- Limit[['min']]
}

# cat(MAX, MIN, "\n")



if(as.character(opt["bin_width"]) != "NA"){
  bin_width <- as.numeric(as.character(opt["bin_width"]))
}else if(as.character(opt["bin_width"]) == "round"){
  bin_width <- round((MAX - MIN) / 30)
  if(bin_width == 0){
    bin_width <- 1
  }
}else{
  bin_width <- (MAX - MIN) / 30
}
b <- seq(Limit[['min']]-bin_width, Limit[['max']]+bin_width, by=bin_width)
l <- lapply(DATA, hist, breaks = b, plot=FALSE)
mids <- unique(unlist(lapply(l, function(x)x$mids)))
densities <- lapply(l, function(x)x$density[match(x=mids, table=x$mids, nomatch=NA)]);
breaks <- unique(unlist(lapply(l, function(x)x$breaks)))
h <- do.call(rbind, densities)
colnames(h) <- format(head(breaks, -1)+bin_width, digits = 3)


num <- as.numeric(colnames(h))
index_col <- colnames(h)[num  >= MIN & num  <= MAX]
h <- h[,index_col]


pdf(NULL)
par(oma=c(0,0,0,0), mar=c(4,4,2,1))
png(filename=FILE_out, width=14, height=11,  units="cm", res = 200)
barplot(h, beside=TRUE, col=Colors, xlab=as.character(opt["xlab"]), ylab="Probability", main=as.character(opt["title"]))
legend(LEGEND_LOC, legend=NAMEs, lwd=5, seg.len = 1.5,col=Colors, bty='n')
if(length(FILE_ins) == 2){
  if(NUM_sampling != "NA"){
    # sampling methodの場合は最初の200samplingのデータだけを用いる
    pval <- (wilcox.test(DATA[[1]][1:200], DATA[[2]][1:200], alternative = "two.sided"))$p.value
  }else{
    pval <- (wilcox.test(DATA[[1]], DATA[[2]], alternative = "two.sided"))$p.value
  }
  
  ratio <- median(DATA[[1]], na.rm = TRUE) / median(DATA[[2]], na.rm = TRUE)

  title(sub=paste0("Pval = ",  format(pval, digits = 5), "  Ratio of average = ", format(ratio, digits = 3)), bty='n', cex.sub=0.8)
}
dummy <- dev.off()


