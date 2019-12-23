#!/usr/bin/Rscript
# Drawing multiple HiC contactmaps very quickly
# output png : [name].png 
# output matrix : [name].matrix

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="directory of matrices located"),
  make_option(c("-o", "--out"),help="output directory"),
  make_option(c("--location"), default="NA", help="file describe target loci"),
  make_option(c("--distance"), default=FALSE, help="normalized by distance curve"),
  make_option(c("--normalize"), default="NA", help="NA, average: average will be 1, probability: score were divided by total read"),
  make_option(c("--moving_average"), default=0, help="number of merging bin for moving average calculation"),
  make_option(c("--na"), default="NA", help="how to treat na value. min, na, ave, zero. min replace with minimum value. ave take average of same distance, zero replace to zero"),
  make_option(c("--zero"), default="NA", help="how to treat 0 value. min, na, ave. min replace with minimum value. ave take average of same distance"),
  make_option(c("--matrix"), default="NULL", help="directory for output matrix. NULL for not output"),
  make_option(c("--color"), default="matlab", help="color matlab or gentle, blue or red"),
  make_option(c("--unit"), default="v", help="unit to define score threshold p:percent or v:value"),
  make_option(c("--min"), default="NULL", help="minimum score for drawing"),
  make_option(c("--max"), default="NULL", help="maximu score for drawing"),
  make_option(c("--width"), default="NULL", help="width of output figure"),
  make_option(c("--height"), default="NULL", help="height of output figure"),
  make_option(c("--linerColor"), default=FALSE, help="use linear color scale"),
  make_option(c("--triangle"), default="FALSE", help="plot only half of triangle")
)
opt <- parse_args(OptionParser(option_list=option_list))

pallete <- c("#00007F", "blue", "#007FFF", 
             "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
if(as.character(opt["color"]) == "gentle"){
  suppressWarnings(suppressMessages(library("RColorBrewer")))
  pallete <- rev(brewer.pal(10, "RdBu"))  # gentle color
}
if(as.character(opt["color"]) == "red"){
  suppressWarnings(library("RColorBrewer"))
  pallete <- brewer.pal(9, "Reds")
}
if(as.character(opt["color"]) == "blue"){
  suppressWarnings(library("RColorBrewer"))
  pallete <- brewer.pal(9, "Blues")
}
if(as.character(opt["color"]) == "yellow"){
  pallete <- rev(c(rgb(0,0,0), rgb(0.13604190193,0.06025908266,0.0358534453625), rgb(0.236596675144,0.0854393225134,0.0663292063043), 
                   rgb(0.345308379611,0.104251700818,0.0859523214123), rgb(0.459650765971,0.118262288197,0.10370568771), 
                   rgb(0.57886917165,0.12712196061,0.121683026697), rgb(0.694333453232,0.147366583596,0.135554156099), 
                   rgb(0.757966046005,0.254807778863,0.115414575242), rgb(0.821231461948,0.346641504267,0.0789177160322), 
                   rgb(0.87357362974,0.441355114094,0.0184645577569), rgb(0.887325704879,0.554597071914,0.013882477102), 
                   rgb(0.895705175617,0.661630633721,0.0133922831093), rgb(0.898208608643,0.76566848057,0.0180193216995), 
                   rgb(0.923430940957,0.85551584532,0.291204862021), rgb(0.977010668325,0.92623393319,0.660287806118), rgb(1,1,1)),space="Lab")
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
D_location <- fread(FILE_location)
if(!("chr2" %in% colnames(D_location))){
  D_location[,"chr2"] = D_location[,"chr1"]
  D_location[,"start2"] = D_location[,"start1"]
  D_location[,"end2"] = D_location[,"end1"]
}

checkDIRpath <- function(DIR){
  if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
    DIR <- paste0(DIR, "/")
  }
  DIR
}
DIR_in <- checkDIRpath(as.character(opt["in"]))
DIR_out <- checkDIRpath(as.character(opt["out"]))
DIR_matrix <- checkDIRpath(as.character(opt["matrix"]))


for(cc in D_location %>% distinct(chr1) %>% pull(chr1) %>% as.character()){
  map <- readRDS(paste0(DIR_in, cc, ".rds"))
  map <- ifelse(is.infinite(map), NA, map)
  r <- rownames(map)
  LocList <- strsplit(r, ":")
  LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
  
  # Normalization
  Normalization = as.character(opt["normalize"])
  if(Normalization == "average"){
    d <- nrow(map)
    map <- map / sum(map, na.rm=TRUE) * d * d
  }else if(Normalization == "probability"){
    map <- map / sum(map, na.rm=TRUE)
  }
  
  if(eval(parse(text=as.character(opt["distance"])))){
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
    map <- ifelse(map_expect == 0, NA, map / map_expect)
    rm(map_expect)
  }
  
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
    map.extract <- map[Region1, Region2]
    
    # moving average
    N_moving_average <- as.numeric(as.character(opt["moving_average"]))
    if(N_moving_average > 0){
      map_new <- map.extract
      for (i in (1:(nrow(map.extract)))){
        for (j in (1:(ncol(map.extract)))){
          map_new[i,j] = mean(map.extract[max(1,i-N_moving_average):min(nrow(map.extract),i+N_moving_average),
                                          max(1,j-N_moving_average):min(ncol(map.extract),j+N_moving_average)], na.rm = TRUE)
        }
      }
      map.extract <- map_new
      rm(map_new)
    }
    
    if(as.character(opt["unit"]) == "p"){
      NUM <- sort(as.numeric(map.extract), na.last = NA)
      if(as.character(opt["min"]) == "NULL"){
        Min <- min(map.extract, na.rm=TRUE)
      }else{
        pct <- as.numeric(as.character(opt["min"]))
        if(pct > 1){
          pct <- pct / 100
        }
        rank <- round(length(NUM)*pct)
        if(rank == 0){
          rank = rank +1
        }
        Min <- NUM[rank]
      }
      if(as.character(opt["max"]) == "NULL"){
        Max <- max(map.extract, na.rm=TRUE)
      }else{
        pct <- as.numeric(as.character(opt["max"]))
        if(pct > 1){
          pct <- pct / 100
        }
        rank <- round(length(NUM)*pct)
        if(rank == 0){
          rank = rank +1
        }
        Max <- NUM[rank]
      }
    }else if(as.character(opt["unit"]) == "v"){
      if(as.character(opt["min"]) == "NULL"){
        Min <- min(map.extract, na.rm=TRUE)
      }else{
        Min <- as.numeric(as.character(opt["min"]))
      }
      if(as.character(opt["max"]) == "NULL"){
        Max <- max(map.extract, na.rm=TRUE)
      }else{
        Max <- as.numeric(as.character(opt["max"]))
      }
    }
    
    
    # replace Na
    if(as.character(opt["na"]) == "min"){
      MinValue <- min(map.extract, na.rm=TRUE)
      map.extract <- ifelse(is.na(map.extract), MinValue, map.extract)
    }else if(as.character(opt["na"]) == "zero"){
      map.extract <- ifelse(is.na(map.extract), 0, map.extract)
    }else if(as.character(opt["na"]) == "ave"){
      index <- which(is.na(map.extract), arr.ind = TRUE)
      if(length(index) > 0){
        index <- index[index[,1] != 1 & index[,1] != nrow(map.extract) & index[,2] != 1 & index[,2] != ncol(map.extract),]
        estimateNa <- function(x, y){
          v <- mean(c(map.extract[x-1,y+1], map.extract[x+1,y-1]))
          if(is.na(v)){
            v <- mean(c(map.extract[x-1,y], map.extract[x+1,y]))
          }
          if(is.na(v)){
            v <- mean(c(map.extract[x,y+1], map.extract[x,y-1]))
          }
          map.extract[x,y] <<- v
        }
        dummy <- mapply(estimateNa, as.integer(index[,1]), as.integer(index[,2]))
      }
    }
    
    # replace 0
    if(as.character(opt["zero"]) == "min"){
      tmp <- map.extract[map.extract != 0]
      MinValue <- min(tmp, na.rm=TRUE)
      map.extract <- ifelse(map.extract ==0, MinValue, map.extract)
    }else if(as.character(opt["zero"]) == "ave"){
      index <- which(map==0, arr.ind = TRUE)
      if(length(index) > 0){
        index <- index[index[,1] > 1 & index[,1] < nrow(map.extract) & index[,2] > 1 & index[,2] < ncol(map.extract),]
        estimateZero <- function(x, y){
          map.extract[x,y] <<- mean(c(map.extract[x-1,y+1], map.extract[x+1,y-1]))
        }
        dummy <- mapply(estimateZero, as.integer(index[,1]), as.integer(index[,2]))
      }
    }
    
    if(DIR_matrix != "NULL/"){
      write.table(map.extract, file=paste0(DIR_matrix, NAME, ".matrix"), quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA)
    }
    
    map.conv <- TakeMiddleV(map.extract, Min, Max)
    
    ### draw only half
    if(eval(parse(text=as.character(opt["triangle"])))){
      half <- lower.tri(map.conv, diag=TRUE)
      map.conv[half] <- NA
    }
    
    
    if(as.character(opt["width"]) == "NULL"){
      width <- nrow(map)
    }else{
      width <- as.numeric(as.character(opt["width"]))
    }
    if(as.character(opt["height"]) == "NULL"){
      height <- width / ncol(map.conv) * nrow(map.conv)
    }else{
      height <- as.numeric(as.character(opt["height"]))
    }
    
    # Colorの調整
    if(eval(parse(text=opt["linerColor"]))){
      bk <- seq(Min, Max, length.out=100)
    }else{
      # t <- (Max - Min) * 0.7 + Min
      # bk <- unique(c(seq(Min, t, length.out=80), lseq(t, Max, length.out = 30), Max+1))
      tmp <- as.numeric(map.conv)
      tmp <- tmp[!is.na(tmp)]
      T95 <- sort(tmp)[round(length(tmp)*0.95)]
      bk <- unique(round(c(seq(Min, T95, length.out=95), lseq(T95, Max+1, length.out=5)), digits=2))
    }
    
    # t <- (Max - Min) * 0.8 + Min
    # bk <- unique(round(c(seq(Min, t, length.out=90), lseq(t, Max+1, length.out = 10)), digits = 2))
    map.cat <- matrix(as.integer(cut(map.conv, breaks = bk, include.lowest = TRUE)), nrow = nrow(map.conv))
    colors <- colorRampPalette(pallete)(length(bk))
    colors <- colors[min(map.cat, na.rm=TRUE):max(map.cat, na.rm=TRUE)]
    
    FILE_OUT <- paste0(DIR_out, NAME, ".png")
    png(file=FILE_OUT, width=width, height=height, units="px", bg="white")
    par(oma=c(0,0,0,0), mar=c(0,0,0,0))
    image(Transform(map.cat), col=colors, axes=F)
    dummy <- dev.off()
  }
  
  dummy <- sapply(1:nrow(D_table), DrawMap)
}

