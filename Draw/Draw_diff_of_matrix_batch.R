#!/usr/bin/Rscript
# Drawing multiple HiC contactmaps very quickly
# output png : [name].png 
# output matrix : [name].matrix

suppressWarnings(suppressMessages(library("spatstat")))
suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-a", "--dir1"),help="directory of matrix 1"),
  make_option(c("-b", "--dir2"),help="directory of matrix 2"),
  make_option(c("-o", "--out"),help="output directory"),
  make_option(c("--location"), default="NA", help="file describe target loci. excel file is okay"),
  make_option(c("--draw_area"), default="NA", help="instead of drawing input file area, draw surrounded xxx area (upstream XXX, downstream XXX) from center"),
  make_option(c("--sufix"), default="NA", help="sufix for output file name."),
  make_option(c("--normalize"), default="NA", help="NA, average: average will be 1, probability: score were divided by total read"),
  make_option(c("--adjust"), default=TRUE, help="adjust total read of file2 to file1"),
  make_option(c("--extract_first"), default=FALSE, help="In default, entire chromosome were considered for adjusting. With this option, Extract target area first then normalize"),
  make_option(c("--blur"), default=TRUE, help="output png with blur"),
  make_option(c("--sigma"), default=".5", help="sigma value of Blur"),
  make_option(c("--matrix"), default="NULL", help="directory for output matrix. NULL for not output"),
  make_option(c("--color"), default="matlab", help="color matlab or gentle, blue or red"),
  make_option(c("--unit"), default="v", help="unit to define score threshold p:percent or v:value"),
  make_option(c("-t", "--threshold"), default="NULL", help="min and max threshold"),
  make_option(c("--method"), default="diff", help="calculation method. diff or ratio"),
  make_option(c("--width"), default="NULL", help="width of output figure"),
  make_option(c("--height"), default="NULL", help="height of output figure"),
  make_option(c("--linerColor"), default=FALSE, help="use linear color scale"),
  make_option(c("--triangle"), default="FALSE", help="plot only half of triangle"),
  make_option(c("--circle"), default="NULL", help="location pairs to draw circles on output")
)
opt <- parse_args(OptionParser(option_list=option_list))


# 色分け
suppressWarnings(suppressMessages(library("RColorBrewer")))
if(as.character(opt["color"]) == "gentle"){
  pallete <- rev(brewer.pal(10, "RdBu"))  # gentle color
}
if(as.character(opt["color"]) == "blight"){
  pallete <- c(rev(brewer.pal(9, "Blues"))[-9], brewer.pal(9, "Reds"))
}


suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))

TakeMiddleV <- function(mat, minVal, maxVal){
  mat.new <- ifelse(mat < minVal, minVal, mat)
  ifelse(mat.new > maxVal, maxVal, mat.new)
}


Transform <- function(mat){
  d = dim(mat)[1]
  n = dim(mat)[2]
  mat.rev = t(mat[d+1-c(1:d), ])
  mat.rev
}

lseq <- function(from=1, to=100000, length.out=6) {
  if(from == to){
    from
  }else{
    exp(seq(log(from), log(to), length.out = length.out))
  }
}

### location file
# chr1, start1, end1, chr2, start2, end2, name
FILE_location <- as.character(opt["location"])
if(substring(FILE_location, nchar(FILE_location)-3, nchar(FILE_location)) == "xlsx"){
  library("xlsx")
  D_location <- read.xlsx(FILE_location, sheetIndex=1, header=TRUE, stringsAsFactors = FALSE)
}else{
  D_location <- fread(FILE_location, header=TRUE, sep="\t", stringsAsFactors = FALSE) %>% as.data.frame()
}


if(!"chr1" %in% colnames(D_location)){
  cat("chr1 should be specified")
  q()
}
if(!"start1" %in% colnames(D_location)){
  cat("start1 should be specified")
  q()
}
if(!"end1" %in% colnames(D_location)){
  cat("end1 should be specified")
  q()
}

### もしnameがなかったら上から順にidを付ける
if(!"name" %in% colnames(D_location)){
  D_location <- D_location %>% mutate(name=row_number())
}

### name sufix
NAME_SUFIX <- as.character(opt["sufix"])
if(NAME_SUFIX != "NA"){
  D_location$name = paste0(D_location$name, NAME_SUFIX)
}

if(!("chr2" %in% colnames(D_location))){
  D_location[,"chr2"] = D_location[,"chr1"]
  D_location[,"start2"] = D_location[,"start1"]
  D_location[,"end2"] = D_location[,"end1"]
}

RANGE_DRAW <- as.character(opt["draw_area"])
if(RANGE_DRAW != "NA"){
  RANGE_DRAW <- as.numeric(RANGE_DRAW)
  D_location <- D_location %>% mutate(center1=as.integer((start1 + end1)/2), center2=as.integer((start2 + end2)/2))
  D_location <- D_location %>% mutate(start1=center1 - RANGE_DRAW, end1=center1 + RANGE_DRAW, start2 = center2 - RANGE_DRAW, end2 = center2 + RANGE_DRAW)
}else{
  draw_range_min <- min(D_location %>% mutate(area=end1-start1) %>% pull(area))
  if(draw_range_min < 20000){
    cat("Some of the draw range are too small:", draw_range)
    q()
  }
}



checkDIRpath <- function(DIR){
  if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
    DIR <- paste0(DIR, "/")
  }
  DIR
}
DIR_in1 <- checkDIRpath(as.character(opt["dir1"]))
DIR_in2 <- checkDIRpath(as.character(opt["dir2"]))
DIR_out <- checkDIRpath(as.character(opt["out"]))
DIR_matrix <- checkDIRpath(as.character(opt["matrix"]))



FILE_circle <- as.character(opt["circle"])
if(FILE_circle != "NULL"){
  DATA_circle <- read.table(as.character(opt["circle"]), header=F, sep="\t", check.names = F)
}

for(cc in D_location %>% distinct(chr1) %>% pull(chr1) %>% as.character()){
  map1 <- readRDS(paste0(DIR_in1, cc, ".rds"))
  map1 <- ifelse(is.infinite(map1), NA, map1)
  r1 <- rownames(map1)
  LocList <- strsplit(r1, ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  
  map2 <- readRDS(paste0(DIR_in2, cc, ".rds"))
  map2 <- ifelse(is.infinite(map2), NA, map2)

  
  D_table <- D_location %>% filter(chr1 == cc)
  
  DrawMap <- function(i){
    NAME <- D_table[i, "name"]
    CHR1 <- D_table[i, "chr1"]
    START1 <- D_table[i, "start1"]
    END1 <- D_table[i, "end1"]
    Region1 <- which((as.character(LocMatrix[,1]) == CHR1) & (as.numeric(LocMatrix[,3]) >= START1) & (as.numeric(LocMatrix[,2]) <= END1))
    CHR2 <- D_table[i, "chr2"]
    START2 <- D_table[i, "start2"]
    END2 <- D_table[i, "end2"]
    Region2 <- which((as.character(LocMatrix[,1]) == CHR2) & (as.numeric(LocMatrix[,3]) >= START2) & (as.numeric(LocMatrix[,2]) <= END2))
    
    map1_new <- map1
    map2_new <- map2
    
    ### 先に領域を絞る
    if(eval(parse(text=opt["extract_first"]))){
      map1_new.tmp <- map1_new[Region1, Region2]
      r2 <- rownames(map2_new)
      Region1 <- intersect(rownames(map1_new.tmp), r2)
      Region2 <- intersect(colnames(map1_new.tmp), r2)
      map2_new <- map2_new[Region1, Region2]
      map1_new <- map1_new[Region1, Region2]
      rm(map1_new.tmp)
    }
    
    # map2の合計値をmap1の合計値と同じになるように調整する
    if(eval(parse(text=opt["adjust"]))){
      TotalRead <- sum(map1_new, na.rm=TRUE)
      map2_new <- map2_new / sum(map2_new, na.rm=TRUE) * TotalRead
    }
    
    # map1とmap2それぞれ、平均値を１にする
    if(eval(parse(text=opt["normalize"]))){
      d <- nrow(map1)
      map1_new <- map1_new / sum(map1_new, na.rm=TRUE) * d * d
      map2_new <- map2_new / sum(map2_new, na.rm=TRUE) * d * d
    }
    
    
    # 領域を絞る（通常の順序)
    if(!eval(parse(text=opt["extract_first"]))){
      map1_new.tmp <- map1_new[Region1, Region2]
      r2 <- rownames(map2_new)
      Region1 <- intersect(rownames(map1_new.tmp), r2)
      Region2 <- intersect(colnames(map1_new.tmp), r2)
      map2_new <- map2_new[Region1, Region2]
      map1_new <- map1_new[Region1, Region2]
      rm(map1_new.tmp)
    }
    
    
    #=========================================================
    # 2つのmapの違いを計算する
    #=========================================================
    if(as.character(opt["method"]) == "diff"){
      mat.diff <- map1_new - map2_new
    }else if(as.character(opt["method"]) == "ratio"){
      mat.diff <- log2(map1_new / map2_new)
      mat.diff <- ifelse(is.infinite(mat.diff), NA, mat.diff)
    }
    rm(map1_new, map2_new)
    
    if(DIR_matrix != "NULL/"){
      write.table(mat.diff, file=paste0(DIR_matrix, NAME, ".matrix"), quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA)
    }
    
    #=========================================================
    # blur image
    #=========================================================
    if(eval(parse(text=opt["blur"]))){
      sigma <- as.numeric(as.character(opt["sigma"]))
      t <- blur(as.im(mat.diff), sigma=sigma, bleed=FALSE)
      mat.diff <- t$v
      rownames(mat.diff) <- r1[Region1]
      colnames(mat.diff) <- r1[Region2]
    }
    

    #=========================================================
    # 閾値で値を区切る
    #=========================================================
    if(as.character(opt["unit"]) == "p"){
      NUM <- sort(abs(as.numeric(mat.diff)), na.last = NA)
      if(as.character(opt["threshold"]) == "NULL"){
        Threshold <- max(NUM, na.rm=TRUE)
      }else{
        pct <- as.numeric(as.character(opt["threshold"]))
        if(pct > 1){
          pct <- pct / 100
        }
        rank <- round(length(NUM)*pct)
        if(rank == 0){
          rank = rank +1
        }
        Threshold <- NUM[rank]
      }
    }else if(as.character(opt["unit"]) == "v"){
      if(as.character(opt["threshold"]) == "NULL"){
        Threshold <- max(abs(mat.diff), na.rm=TRUE)
      }else{
        Threshold <- as.numeric(as.character(opt["threshold"]))
      }
    }
    
    mat.diff <- ifelse(mat.diff < (Threshold * -1), Threshold * -1, mat.diff)
    mat.diff <- ifelse(mat.diff > Threshold, Threshold, mat.diff)
    
  
    tmp <- as.numeric(mat.diff)
    tmp <- tmp[!is.na(tmp)]
    bk <- seq(-Threshold, Threshold, length.out=100)
    mat.diff <- matrix(as.integer(cut(mat.diff, breaks = bk, include.lowest = TRUE)), nrow = nrow(mat.diff), ncol=ncol(mat.diff))
    
    colors <- colorRampPalette(pallete)(length(bk))
    colors <- colors[min(mat.diff, na.rm=TRUE):max(mat.diff, na.rm=TRUE)]
    
    
    ### draw only half
    if(eval(parse(text=as.character(opt["triangle"])))){
      half <- lower.tri(mat.diff, diag=TRUE)
      mat.diff[half] <- NA
    }
    
    
    if(as.character(opt["width"]) == "NULL"){
      width <- nrow(mat.diff)
    }else{
      width <- as.numeric(as.character(opt["width"]))
    }
    if(as.character(opt["height"]) == "NULL"){
      height <- width / ncol(mat.diff) * nrow(mat.diff)
    }else{
      height <- as.numeric(as.character(opt["height"]))
    }
    
    
    FILE_OUT <- paste0(DIR_out, NAME, ".png")
    png(file=FILE_OUT, width=width, height=height, units="px", bg="white")
    par(oma=c(0,0,0,0), mar=c(0,0,0,0))
    image(Transform(mat.diff), col=colors, axes=F)
    
    ### Circle
    if(FILE_circle != "NULL"){
      for(i in 1:nrow(DATA_circle)){
        if(DATA_circle[i,1] %in% Region1 & DATA_circle[i,2] %in% Region2){
          par(new=T)
          plot(which(DATA_circle[i,1] == Region1), nrow(mat.diff) - which(DATA_circle[i,2] == Region2)+1, pch=21, xlim=c(1,ncol(mat.diff)), 
               ylim=c(1,nrow(mat.diff)), xaxs="i", yaxs="i", cex=2, axes=F, col='black', lwd=2.5)
        }
      }
    }
    dummy <- dev.off()
  }
  
  dummy <- sapply(1:nrow(D_table), DrawMap)
}

