#!/usr/bin/Rscript
# Draw arch of associations

#=============================================
# Input file should have at least following column header
# start1, end1, start2, end2
#=============================================

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"), default="NA", help="input file"),
  make_option(c("-o", "--out"), default="NA", help="output image file"),
  make_option(c("--chr"), default="NA", help="chromosome name"),
  make_option(c("--start"), default="NA", help="start position for drawing"),
  make_option(c("--end"), default="NA", help="end position for drawing"),
  make_option(c("--min"), default=0, help="minimum distance to drawing"),
  make_option(c("--bottom"), default="TRUE", help="draw bottom (TRUE) or not (FALSE)"),
  make_option(c("--color"), default="grey30", help="colors")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))

FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
CHR <- as.character(opt["chr"])
START <- as.numeric(as.character(opt["start"]))
END <- as.numeric(as.character(opt["end"]))
MIN_dis <- as.numeric(as.character(opt["min"]))
COLOR <- as.character(opt["color"])
COLOR <- adjustcolor(COLOR, alpha.f = 0.4)
FLAG_bottom <- eval(parse(text=as.character(opt["bottom"])))

getArchValues <- function(p1, p2){
  Px <- seq(p1, p2, length.out = 100)
  radius <- abs(p2 - p1)/2
  middle <- (p1 + p2)/2
  Py <- sqrt(radius**2 - (Px-middle)**2)
  if(FLAG_bottom){
    Py <- Py * -1
  }
  out <- data.frame(x=Px, y=Py, stringsAsFactors = FALSE)
  out
}

drawConnecLine <- function(s1,e1, s2, e2){
  sp1 <- getArchValues(s1, e2)
  sp2 <- getArchValues(e1, s2)
  xx <- c(sp1$x, rev(sp2$x))
  yy <- c(sp1$y, rev(sp2$y))
  
  polygon(xx, yy, border = NA, col=COLOR)
}

drawEach <- function(i){
  D_sub <- D_table[i,]
  drawConnecLine(D_sub$start1, D_sub$end1, D_sub$start2, D_sub$end2)
}


D_table <- read.table(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE)
if(sum(D_table %>% colnames() %in% c("start1", "end1", "start2", "end2")) != 4){
  print("ERROR: at least start1, end1, start2, end2 are required for input file!")
  q()
}

index_reverse <- which(D_table$start1 > D_table$start2)
D_table[index_reverse, c("start1", "end1", "start2", "end2")] = D_table[index_reverse, c("start2", "end2", "start1", "end1")]


if(CHR != "NA"){
  D_table <- D_table %>% filter(chr1==CHR & chr2==CHR)
}
if(START != "NA"){
  D_table <- D_table %>% filter(start1 > START)
}else{
  START = D_table %>% pull(start1) %>% min()
}
if(END != "NA"){
  D_table <- D_table %>% filter(end2 < END)
}else{
  END = D_table %>% pull(end2) %>% max()
}
D_table <- D_table %>% filter(start2-end1 > MIN_dis)

if(nrow(D_table) > 2000){
  print("ERROR: too much lines to draw!")
  q()
}

png(file=FILE_out, width=800, height=400, units="px", bg="transparent")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
Ymax <- (END-START)/2
if(FLAG_bottom){
  Ymax <- Ymax * -1
}
plot(c(START, END), c(0, Ymax), type='n', xlab="", ylab="", xaxs="i", yaxs="i",bty="n", axes=F)
dummy <- sapply(1:nrow(D_table), drawEach)
dummy <- dev.off()



