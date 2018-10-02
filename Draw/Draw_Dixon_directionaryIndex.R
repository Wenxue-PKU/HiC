# DI indexを計算する

# .libPaths(.libPaths()[c(2,3,1)])
suppressPackageStartupMessages(library(HiTC))

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="rds file"),
  make_option(c("-t", "--hitc"),help="rds file of HiTC"),
  make_option(c("-o", "--out"), default="NA", help="output png file"),
  make_option(c("--location"), default="NA", help="output location information as file"),
  make_option(c("-r", "--resolution"), default="40kb", help="resolution of analysis"),
  make_option(c("--window"), default="NA", help="window size to DI calculation"),
  make_option(c("--chr"),help="chromosome name all for all chromosome"),
  make_option(c("--start"), default="1", help="start position"),
  make_option(c("--end"), default="all", help="end position. all for end of the chromosome"),
  make_option(c("--width"), default="NA", help="width of picture"),
  make_option(c("--height"), default="NA", help="height of picture"),
  make_option(c("--reverse"), default="FALSE", help="reverse score plus to minus")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_object <- as.character(opt["in"])
FILE_matrix <- sub(".rds", ".matrix.gz", FILE_object)
HiTC_object <- as.character(opt["hitc"])

# convert matrix file to HiTC format
# get path of program directory
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
ConvetProgram <- paste(sep="/", script.basename, "../TAD/Convert_matrix2HiTC.pl")

FILE_matrix2 <- paste(HiTC_object, "_matrix_HiTC.txt", sep="")


if(file.exists(HiTC_object)){
  hic <- readRDS(HiTC_object)
}else{
  if(eval(parse(text=as.character(opt["reverse"])))){
    conversionCommand <- paste("perl", ConvetProgram, "-i", FILE_matrix, "-o", FILE_matrix2, "-p")
  }else{
    conversionCommand <- paste("perl", ConvetProgram, "-i", FILE_matrix, "-o", FILE_matrix2)
  }
  system(conversionCommand)
  hic <- import.my5C(FILE_matrix2)
  saveRDS(hic, HiTC_object)
  file.remove(FILE_matrix2)
}


map <- readRDS(FILE_object)

r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
if(as.character(opt["end"]) == "all"){
  SameChromosome <- which(as.character(LocMatrix[,1]) == CHR)
  END <- max(as.numeric(LocMatrix[SameChromosome,3]));
}else{
  END <- as.numeric(as.character(opt["end"]))
}


if(as.character(opt["window"]) != "NA"){
  window <- as.numeric(as.character(opt["window"]))
}else{
  window <- switch(as.character(opt["resolution"]), "5kb"=40000, "10kb"=60000, "20kb"=1000000, "40kb"=2000000)
}

if(CHR == "I" || CHR=="chr1"){DI <- directionalityIndex(hic$chr1chr1, winup = window, windown = window)
} else if(CHR == "II" || CHR=="chr2"){DI <- directionalityIndex(hic$chr2chr2, winup = window, windown = window)
} else if(CHR == "III" || CHR=="chr3"){DI <- directionalityIndex(hic$chr3chr3, winup = window, windown = window)
} else if(CHR=="chr4"){ DI <- directionalityIndex(hic$chr4chr4, winup = window, windown = window)
} else if(CHR=="chr5"){ DI <- directionalityIndex(hic$chr5chr5, winup = window, windown = window)
} else if(CHR=="chr6"){ DI <- directionalityIndex(hic$chr6chr6, winup = window, windown = window)
} else if(CHR=="chr7"){ DI <- directionalityIndex(hic$chr7chr7, winup = window, windown = window)
} else if(CHR=="chr8"){ DI <- directionalityIndex(hic$chr8chr8, winup = window, windown = window)
} else if(CHR=="chr9"){ DI <- directionalityIndex(hic$chr9chr9, winup = window, windown = window)
} else if(CHR=="chr10"){ DI <- directionalityIndex(hic$chr10chr10, winup = window, windown = window)
} else if(CHR=="chr11"){ DI <- directionalityIndex(hic$chr11chr11, winup = window, windown = window)
} else if(CHR=="chr12"){ DI <- directionalityIndex(hic$chr12chr12, winup = window, windown = window)
} else if(CHR=="chr13"){ DI <- directionalityIndex(hic$chr13chr13, winup = window, windown = window)
} else if(CHR=="chr14"){ DI <- directionalityIndex(hic$chr14chr14, winup = window, windown = window)
} else if(CHR=="chr15"){ DI <- directionalityIndex(hic$chr15chr15, winup = window, windown = window)
} else if(CHR=="chr16"){ DI <- directionalityIndex(hic$chr16chr16, winup = window, windown = window)
} else if(CHR=="chr17"){ DI <- directionalityIndex(hic$chr17chr17, winup = window, windown = window)
} else if(CHR=="chr18"){ DI <- directionalityIndex(hic$chr18chr18, winup = window, windown = window)
} else if(CHR=="chr19"){ DI <- directionalityIndex(hic$chr19chr19, winup = window, windown = window)
} else if(CHR=="chr20"){ DI <- directionalityIndex(hic$chr20chr20, winup = window, windown = window)
} else if(CHR=="chr21"){ DI <- directionalityIndex(hic$chr21chr21, winup = window, windown = window)
} else if(CHR=="chr22"){ DI <- directionalityIndex(hic$chr22chr22, winup = window, windown = window)
} else if(CHR=="chrX"){ DI <- directionalityIndex(hic$chrXchrX, winup = window, windown = window)
} else if(CHR=="chrY"){ DI <- directionalityIndex(hic$chrYchrY, winup = window, windown = window)
} else if(CHR=="chr23"){ DI <- directionalityIndex(hic$chr23chr23, winup = window, windown = window)
}

# 指定された場所
SameChromosomeLoc <- LocMatrix[which(as.character(LocMatrix[,1]) == CHR),]
Region <- which((as.numeric(SameChromosomeLoc[,3]) >= START) & (as.numeric(SameChromosomeLoc[,2]) <= END))
DI <- DI[Region]

# boundaryを計算
boundaries = c()
for(i in 2:length(DI)){
  if(DI[i-1] < 0 & DI[i] > 0){
    boundaries <- c(boundaries, i-1)
  }
}

FILE_OUT <- as.character(opt["out"])
if(FILE_OUT != "NA"){
  if(as.character(opt["height"])=="NA"){
    if(sum((grep("\\.eps$", FILE_OUT))) == 1){
      h <- 0.8
    }else{
      h <- 50
    }
  }else{
    h <- as.numeric(as.character(opt["height"]))
  }
  if(as.character(opt["width"])=="NA"){
    if(sum((grep("\\.eps$", FILE_OUT))) == 1){
      w <- 4
    }else{
      w <- 1000
    }
  }else{
    w <- as.numeric(as.character(opt["width"]))
  }
  if(sum((grep("\\.eps$", FILE_OUT))) == 1){
    postscript(file=FILE_OUT, horizontal=FALSE, onefile=FALSE, paper="special", height=h, width=w, family="Helvetica")
  }
  if(sum((grep("\\.png$", FILE_OUT))) == 1){
    png(filename=FILE_OUT, width=w, height=h)
  }
  par(oma=c(0,0,0,0), mar=c(1,0,1,0))
  barplot(DI, col=ifelse(DI>0,"darkred", "darkgreen"), space = 0, xaxs="i", yaxs="i", axes=F, border=NA)
  abline(v=boundaries, col='blue')
  dummy <- dev.off()
}


### 座標を出力する
FILE_location <- as.character(opt["location"])
if(FILE_location != "NA"){
  bd <- rep(0, length(DI))
  bd[boundaries] = 1
  Region <- which((as.character(LocMatrix[,1]) == CHR) & (as.numeric(LocMatrix[,3]) >= START) & (as.numeric(LocMatrix[,2]) <= END))
  OUTPUT <- cbind(LocMatrix[Region,], DI, bd)
  write.table(OUTPUT, file=FILE_location, quote=FALSE, sep="\t", eol="\n", col.names=FALSE, row.names=FALSE, append=FALSE)
}



