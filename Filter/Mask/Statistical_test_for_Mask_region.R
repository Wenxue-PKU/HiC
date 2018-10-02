# Statistics at Mask region

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("--file1"),help="matrix file1"),
  make_option(c("--file2"),help="matrix file2"),
  make_option(c("--bio1"),help="matrix file biological replica1"),
  make_option(c("--bio2"),help="matrix file biological replica2"),
  make_option(c("--adjust_control"), default=FALSE, help="matrix file biological replica2"),
  make_option(c("--FDR"), default="0.05", help="FDR threshold"),
  make_option(c("--mask"), default="NA", help="mask .rds file"),
  make_option(c("--xmin"), default="NA", help="xmin"),
  make_option(c("--xmax"), default="NA", help="xmax"),
  make_option(c("--ymin"), default="NA", help="ymin"),
  make_option(c("--ymax"), default="NA", help="ymax"),
  make_option(c("--title"), default="", help="title of comparison"),
  make_option(c("--image"), default="NA", help="file name of scatter plot"),
  make_option(c("--color_set"), default="RGB", help="RYB: Red Yellow Blue, RGB : Red gray blue"),
  make_option(c("--fill"), default="FALSE", help="fill the image"),
  make_option(c("--pie"), default="NA", help="file name of pie chart"),
  make_option(c("--box"), default="NA", help="file name of box plot"),
  make_option(c("--box_data"), default="NA", help="file name of box plot data file"),
  make_option(c("--distance_distribution"), default="NA", help="file name of distance distribution box plot")
)
opt <- parse_args(OptionParser(option_list=option_list))

# test data
FILE_object1 <- "W:/Data/2017-03-28_HiC_pombe_revising_paper/time8_r1/10kb/ICN/ALL.rds"
FILE_object2 <- "W:/Data/2017-03-28_HiC_pombe_revising_paper/time3_r1/10kb/ICN/ALL.rds"
FILE_object3 <- "W:/Data/2017-03-28_HiC_pombe_revising_paper/time8_r1/10kb/ICN/ALL.rds"
FILE_object4 <- "W:/Data/2017-03-28_HiC_pombe_revising_paper/time8_r2/10kb/ICN/ALL.rds"
FILE_mask <- "W:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-04-12_domain_statistics/time3_r1/Mask_Big_ALL.rds"
TITLE=""
FILE_image <- "W:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/img/time3_time8.png"
FILE_pie <- "W:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/img/pie.png"

FILE_box <- as.character(opt["box"])
FILE_box_data <- as.character(opt["box_data"])
FILE_pie <- as.character(opt["pie"])
FILE_image <- as.character(opt["image"])
FILE_distanceDist <- as.character(opt["distance_distribution"])
TITLE <- as.character(opt["title"])
FILE_mask <- as.character(opt["mask"])
FILE_matrix1 <- as.character(opt["file1"])
FILE_matrix2 <- as.character(opt["file2"])
FILE_matrix3 <- as.character(opt["bio1"])
FILE_matrix4 <- as.character(opt["bio2"])
FILE_object1 <- sub("matrix", "rds", FILE_matrix1)
FILE_object2 <- sub("matrix", "rds", FILE_matrix2)
FILE_object3 <- sub("matrix", "rds", FILE_matrix3)
FILE_object4 <- sub("matrix", "rds", FILE_matrix4)

map1 <- readRDS(FILE_object1)
map2 <- readRDS(FILE_object2)
map_bio1 <- readRDS(FILE_object3)
map_bio2 <- readRDS(FILE_object4)

if(FILE_mask != "NA"){
  mask <- readRDS(FILE_mask)
}else{
  mask <- ifelse(map1 > 0, TRUE, FALSE)
}

r1 <- rownames(map1)
r2 <- rownames(map2)
r3 <- rownames(map_bio1)
r4 <- rownames(map_bio2)
r5 <- rownames(mask)
r.common <- intersect(r1, r2)
r.common <- intersect(r.common, r5)
r.common2 <- intersect(r3, r4)

LocList <- strsplit(r.common, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1

map1 <- map1[r.common, r.common]
map2 <- map2[r.common, r.common]
map_bio1 <- map_bio1[r.common2, r.common2]
map_bio2 <- map_bio2[r.common2, r.common2]
mask <- mask[r.common, r.common]

d <- nrow(map1)
map1 <- map1 / sum(map1, na.rm=TRUE) * d * d
map2 <- map2 / sum(map2, na.rm=TRUE) * d * d
d <- nrow(map_bio1)
map_bio1 <- map_bio1 / sum(map_bio1, na.rm=TRUE) * d * d
map_bio2 <- map_bio2 / sum(map_bio2, na.rm=TRUE) * d * d


target.ratio <- log2(map1/map2)
target.ratio <- ifelse(is.infinite(target.ratio), NA, target.ratio)
target.ave <- (map1 + map2)/2

control.ratio <- as.numeric(log2(map_bio1/map_bio2))
control.ave <- as.numeric((map_bio1 + map_bio2)/2)
control.strange <- is.na(control.ratio) | is.infinite(control.ratio)
control.ratio <- control.ratio[!control.strange]
control.ave <- control.ave[!control.strange]


if(eval(parse(text=as.character(opt["adjust_control"])))){
  index <- sample(length(control.ratio), length(control.ratio)/2, replace = FALSE)
  control.ratio[index] <- control.ratio[index] * -1
}

# # controlのaverageとratioのグラフ
# FILE_control_graph <- "W:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/img/Ave_log2_replica2.png"
# png(filename=FILE_control_graph, width=300, height=300, units="px", bg="white")
# # par(mar=c(4,4,1,1), oma=c(0,0,0,0))
# par(mar=c(0,0,0,0), oma=c(0,0,0,0))
# plot(control.ave, control.ratio, xlab="", ylab="", log="x", main="", xaxs="i", yaxs="i",　pch=20, cex.lab = 1.5, cex.axis = 1.2, cex=0.5,　bty = "n",  ylim=c(-5,5), xlim=c(0.01, 500), axes=F)
# dev.off()


lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

# controlが最低でも1000あるようにカテゴリー分け
control.ave.sort <- sort(control.ave, na.last = NA)
index_max <- length(control.ave.sort)-5000
b <- c(0, lseq(control.ave.sort[5000], control.ave.sort[index_max], 20),1e9)

target.category <- sapply(1:nrow(target.ave), function(i){cut(target.ave[,i], breaks=b, right=FALSE)})
control.category <- cut(control.ave, breaks = b, right=FALSE)

NUM_sampling <- 10000
background <- list()
category_list <- levels(factor(control.category))
for(ca in category_list){
  background[[ca]] <- sample(control.ratio[which(control.category == ca)], NUM_sampling, replace = TRUE)
}

pvalCalc <- function(cate, val){
  if(is.na(val)){
    NA
  }else if(exists(as.character(cate), where=background)){
    sum(val > background[[cate]]) / NUM_sampling
  }else{
    NA
  }
}
pval <- sapply(1:nrow(target.ratio), function(i){mapply(pvalCalc, target.category[i,], target.ratio[i,])})
pval <- ifelse(pval > 0.5, 1-pval, pval)
qval <- p.adjust(pval, method = "BH")

FDR_threshold <- as.numeric(as.character(opt["FDR"]))
pval_threshold <- max(as.numeric(pval[qval <= FDR_threshold]), na.rm = TRUE)

up <- pval <= pval_threshold & target.ratio > 0
down <- pval <= pval_threshold & target.ratio < 0

# # 論文用に細分化する
# down1 <- pval <= pval_threshold & target.ratio < 0 & target.ave < 10
# down2 <- pval <= pval_threshold & target.ratio < 0 & target.ave > 10
# up1 <- pval <= pval_threshold & target.ratio > 0 & target.ave < 30
# up2 <- pval <= pval_threshold & target.ratio > 0 & target.ave > 30


NUM_up <- sum(up[mask])
NUM_down <- sum(down[mask])
NUM_other <- length(up[mask]) - NUM_up - NUM_down
slices <- c(NUM_up, NUM_down, NUM_other)
pct <- round(slices/sum(slices)*100)

cat(paste("UP:", NUM_up, " (", pct[1], " %)\n", sep=""))
cat(paste("DOWN:", NUM_down, " (", pct[2], " %)\n", sep=""))
cat(paste("OTHER:", NUM_other, " (", pct[3], " %)\n\n", sep=""))

t <- summary(target.ratio[mask])
for(i in 1:length(t)){
  cat(names(t[i]), "\t", t[i], "\n")
}



if(as.character(opt["xmin"]) != "NA"){
  xmin <- as.numeric(as.character(opt["xmin"]))
}else{
  xmin <- min(target.ave[mask])
}
if(as.character(opt["xmax"]) != "NA"){
  xmax <- as.numeric(as.character(opt["xmax"]))
}else{
  xmax <- max(target.ave[mask])
}

if(as.character(opt["ymin"]) != "NA"){
  ymin <- as.numeric(as.character(opt["ymin"]))
}else{
  ymin <- min(target.ratio[mask])
}
if(as.character(opt["ymax"]) != "NA"){
  ymax <- as.numeric(as.character(opt["ymax"]))
}else{
  ymax <- max(target.ratio[mask])
}

COLOR_setting <- as.character(opt["color_set"])
if(COLOR_setting  == "RYB"){
  COLOR_set <- c("#FBF7C1" ,"#D72027", "#2C7BB6")
}else if(COLOR_setting == "RGB"){
  COLOR_set <- c(rgb(160,160,160,maxColorValue = 255), rgb(177,1,39,maxColorValue = 255), rgb(85,72,193,maxColorValue = 255))
}


if(FILE_pie != "NA"){
  png(file=FILE_pie, width=150, height=150, units="px", bg="white")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  pie(slices, labels=pct, border=NA, angle=0,
      col=COLOR_set, cex=1.5)
  dummy <- dev.off()
}

if(FILE_box != "NA"){
  # png(file=FILE_box, width=80, height=200, units="px", bg="white")
  postscript(file=FILE_box, horizontal=FALSE, onefile=FALSE, paper="special", height=2, width=1, family="Helvetica")
  par(oma=c(0,0,0,0), mar=c(0,0,0,0))
  boxplot(target.ratio[mask], ylim=c(ymin, ymax), axes=F, col='gray')
  dummy <- dev.off()
}

if(FILE_box_data != "NA"){
  write.table(target.ratio[mask], file=FILE_box_data, quote = FALSE, row.names = TRUE, col.names = NA, sep="\t", eol = "\n")
}



if(FILE_image != "NA"){
  colorList <- matrix(rep(COLOR_set[1], length(up)), nrow=nrow(up))
  colorList[up] <- COLOR_set[2]
  colorList[down] <- COLOR_set[3]
  
  png(file=FILE_image, width=300, height=300, units="px", bg="white")
  # png(file=FILE_image, width=200, height=300, units="px", bg="white")
  if(eval(parse(text=as.character(opt["fill"])))){
    par(oma=c(0,0,0,0), mar=c(0,0,0,0))
    plot(target.ave[mask], target.ratio[mask], pch=20, col=colorList[mask], xlim=c(xmin, xmax), ylim=c(ymin, ymax), bty = "n", axes = F,
         log='x', xlab="", ylab="", xaxs="i", yaxs="i", main="", cex=0.5, cex.axis=1.2)
  }else{
    par(oma=c(4,4,1,1), mar=c(0,0,0,0))
    plot(target.ave[mask], target.ratio[mask], pch=20, col=colorList[mask], xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
         log='x', xlab="Average", ylab="Log2 ratio", xaxs="i", yaxs="i", main="", cex=0.5, cex.axis=1.2)
  }
  dummy <- dev.off()
}

if(FILE_distanceDist != "NA"){
  intraChr <- function(x){
    as.character(LocMatrix[x,1]) == as.character(LocMatrix[,1])
  }
  mask.intraChr <- cbind(sapply(1:length(r.common), intraChr))
  
  map_distance <- map1
  
  NUM_LINE <- length(r.common)
  for(d in 0:(NUM_LINE-1)){
    index1 <- 1:(NUM_LINE - d)
    index2 <- index1 + d
    index3 <- cbind(index1, index2)
    index4 <- cbind(index2, index1)
    Distance <- d * Resolution
    map_distance[index3] <- Distance
    map_distance[index4] <- Distance
  }
  map_distance[!mask.intraChr] <- NA
  
  DATA_for_box <- list()
  DATA_for_box[["UP"]] <- as.numeric(map_distance[mask & up])/1000
  DATA_for_box[["DOWN"]] <- as.numeric(map_distance[mask & down])/1000
  # DATA_for_box[["DOWN1"]] <- as.numeric(map_distance[mask & down1])/1000
  # DATA_for_box[["DOWN2"]] <- as.numeric(map_distance[mask & down2])/1000
  # DATA_for_box[["UP1"]] <- as.numeric(map_distance[mask & up1])/1000
  # DATA_for_box[["UP2"]] <- as.numeric(map_distance[mask & up2])/1000
  
  postscript(file=FILE_distanceDist, horizontal=FALSE, onefile=FALSE, paper="special", height=3, width=3, family="Helvetica")
  par(mar=c(4,3,1,1), oma=c(0,0,0,0))
  boxplot(DATA_for_box, ylim=c(1,5000), col='gray', log="y")
  dummy <- dev.off()
  
  # postscript(file=paste(FILE_distanceDist, "_histogram1.eps", sep=""), horizontal=FALSE, onefile=FALSE, paper="special", height=3, width=3, family="Helvetica")
  # b <- seq(0,6000, length.out=60)
  # par(mar=c(4,3,1,1), oma=c(0,0,0,0))
  # hist(DATA_for_box[["UP"]], col=rgb(0.9,0,0), xlim=c(0,5000), breaks = b, main="", xlab="distance")
  # dummy <- dev.off()
  # 
  # postscript(file=paste(FILE_distanceDist, "_histogram2.eps", sep=""), horizontal=FALSE, onefile=FALSE, paper="special", height=3, width=3, family="Helvetica")
  # par(mar=c(4,3,1,1), oma=c(0,0,0,0))
  # hist(DATA_for_box[["DOWN1"]], col=rgb(0,0,0.9), xlim=c(0,5000), breaks = b, main="", xlab="distance")
  # dummy <- dev.off()
  # 
  # postscript(file=paste(FILE_distanceDist, "_histogram3.eps", sep=""), horizontal=FALSE, onefile=FALSE, paper="special", height=3, width=3, family="Helvetica")
  # par(mar=c(4,3,1,1), oma=c(0,0,0,0))
  # hist(DATA_for_box[["DOWN2"]], col=rgb(0,0.9,0), xlim=c(0,5000), breaks = b, main="", xlab="distance")
  # dummy <- dev.off()

  postscript(file=paste(FILE_distanceDist, "_histogram1v2.eps", sep=""), horizontal=FALSE, onefile=FALSE, paper="special", height=3, width=3, family="Helvetica")
  b <- seq(0,6000, length.out=300)
  par(mar=c(4,3,1,1), oma=c(0,0,0,0))
  hist(DATA_for_box[["UP"]], col=rgb(0.9,0,0), xlim=c(0,1000), breaks = b, main="", xlab="distance")
  dummy <- dev.off()
  
  postscript(file=paste(FILE_distanceDist, "_histogram2v2.eps", sep=""), horizontal=FALSE, onefile=FALSE, paper="special", height=3, width=3, family="Helvetica")
  par(mar=c(4,3,1,1), oma=c(0,0,0,0))
  hist(DATA_for_box[["DOWN"]], col=rgb(0,0.9,0), xlim=c(0,1000), breaks = b, main="", xlab="distance")
  dummy <- dev.off()
}

# # # 論文用のサンプルデータ
# colorList <- rep(COLOR_set[1], length(qval))
# colorList[up] <- COLOR_set[3]
# colorList[down] <- COLOR_set[2]
# png(filename=FILE_image, width=300, height=300, units="px", bg="white")
# par(mar=c(0,0,0,0), oma=c(0,0,0,0))
# plot(target.ave, -target.ratio, xlab="", ylab="", log="x", main="", pch=20,
#      cex.lab = 1.5, cex.axis = 1.2, cex=0.5,  ylim=c(-5,5), xlim=c(0.005, 500), axes=F, col=colorList)
# dummy <- dev.off()

