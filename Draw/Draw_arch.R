#!/usr/bin/Rscript
# Draw arch of associations

#=============================================
# Input file should have at least following column header
# start1, end1, start2, end2
# combのカラムがある場合には、"E-P"だけが色で指定したものの対象になる
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
  make_option(c("--color"), default="grey30", help="colors. 'gentle' to use gentle color"),
  make_option(c("--score"), default="NA", help="target column for color definition"),
  make_option(c("--cairo"), default="TRUE", help="use cairo for output")
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
if(COLOR == "gentle"){
  COLOR <- colorRampPalette(c(rgb(85,72,193, max=255), rgb(88,76,196, max=255), rgb(90,79,199, max=255), rgb(92,83,202, max=255), rgb(94,87,205, max=255), rgb(96,90,208, max=255), rgb(98,94,211, max=255), rgb(100,97,214, max=255), rgb(103,101,216, max=255), rgb(105,104,219, max=255), rgb(107,108,221, max=255), rgb(109,111,224, max=255), rgb(111,114,226, max=255), rgb(114,118,229, max=255), rgb(116,121,231, max=255), rgb(118,124,233, max=255), rgb(120,128,235, max=255), rgb(122,131,237, max=255), rgb(125,134,239, max=255), rgb(127,137,240, max=255), rgb(129,140,242, max=255), rgb(131,144,244, max=255), rgb(134,147,245, max=255), rgb(136,150,246, max=255), rgb(138,153,248, max=255), rgb(140,156,249, max=255), rgb(143,158,250, max=255), rgb(145,161,251, max=255), rgb(147,164,252, max=255), rgb(149,167,253, max=255), rgb(152,169,253, max=255), rgb(154,172,254, max=255), rgb(156,175,254, max=255), rgb(159,177,255, max=255), rgb(161,180,255, max=255), rgb(163,182,255, max=255), rgb(165,184,255, max=255), rgb(168,187,255, max=255), rgb(170,189,255, max=255), rgb(172,191,255, max=255), rgb(174,193,255, max=255), rgb(176,195,254, max=255), rgb(179,197,254, max=255), rgb(181,199,253, max=255), rgb(183,201,253, max=255), rgb(185,203,252, max=255), rgb(187,204,251, max=255), rgb(189,206,250, max=255), rgb(192,207,249, max=255), rgb(194,209,248, max=255), rgb(196,210,246, max=255), rgb(198,211,245, max=255), rgb(200,213,244, max=255), rgb(202,214,242, max=255), rgb(204,215,240, max=255), rgb(206,216,239, max=255), rgb(208,217,237, max=255), rgb(209,218,235, max=255), rgb(211,218,233, max=255), rgb(213,219,231, max=255), rgb(215,219,229, max=255), rgb(217,220,227, max=255), rgb(218,220,224, max=255), rgb(220,221,222, max=255), rgb(222,220,219, max=255), rgb(224,219,216, max=255), rgb(225,218,214, max=255), rgb(227,217,211, max=255), rgb(229,216,208, max=255), rgb(230,215,205, max=255), rgb(232,213,202, max=255), rgb(233,212,199, max=255), rgb(234,210,196, max=255), rgb(236,209,193, max=255), rgb(237,207,190, max=255), rgb(238,205,187, max=255), rgb(239,203,183, max=255), rgb(239,201,180, max=255), rgb(240,199,177, max=255), rgb(241,197,174, max=255), rgb(242,195,171, max=255), rgb(242,193,168, max=255), rgb(242,191,165, max=255), rgb(243,188,161, max=255), rgb(243,186,158, max=255), rgb(243,183,155, max=255), rgb(243,181,152, max=255), rgb(243,178,149, max=255), rgb(243,175,146, max=255), rgb(243,173,142, max=255), rgb(243,170,139, max=255), rgb(242,167,136, max=255), rgb(242,164,133, max=255), rgb(241,161,130, max=255), rgb(241,158,127, max=255), rgb(240,155,124, max=255), rgb(239,152,121, max=255), rgb(239,148,118, max=255), rgb(238,145,115, max=255), rgb(237,142,111, max=255), rgb(235,138,109, max=255), rgb(234,135,106, max=255), rgb(233,131,103, max=255), rgb(232,128,100, max=255), rgb(230,124,97, max=255), rgb(229,120,94, max=255), rgb(227,117,91, max=255), rgb(225,113,88, max=255), rgb(224,109,85, max=255), rgb(222,105,83, max=255), rgb(220,101,80, max=255), rgb(218,97,77, max=255), rgb(216,93,75, max=255), rgb(214,89,72, max=255), rgb(212,84,69, max=255), rgb(209,80,67, max=255), rgb(207,76,64, max=255), rgb(205,71,62, max=255), rgb(202,66,59, max=255), rgb(200,61,57, max=255), rgb(197,56,55, max=255), rgb(194,51,52, max=255), rgb(192,45,50, max=255), rgb(189,39,48, max=255), rgb(186,33,46, max=255), rgb(183,25,44, max=255),
                                 rgb(180,15,41, max=255), rgb(177,1,39, max=255)), space="Lab")( 21 )
  COLOR <- adjustcolor(COLOR, alpha.f = 0.5)
}else{
  COLOR <- c(adjustcolor(COLOR, alpha.f = 0.5), adjustcolor("grey50", alpha.f = 0.2))
}
FLAG_bottom <- eval(parse(text=as.character(opt["bottom"])))
FLAG_cairo <- eval(parse(text=as.character(opt["cairo"])))
TARGET_column <- as.character(opt["score"])

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

drawConnecLine <- function(s1,e1, s2, e2, ccc){
  sp1 <- getArchValues(s1, e2)
  sp2 <- getArchValues(e1, s2)
  xx <- c(sp1$x, rev(sp2$x))
  yy <- c(sp1$y, rev(sp2$y))
  polygon(xx, yy, border = NA, col=ccc)
}

drawEach <- function(i){
  D_sub <- D_table[i,]
  drawConnecLine(D_sub$start1, D_sub$end1, D_sub$start2, D_sub$end2, D_sub$z)
}


D_table <- read.table(FILE_in, header=TRUE, sep="\t", stringsAsFactors = FALSE)
if(sum(D_table %>% colnames() %in% c("start1", "end1", "start2", "end2")) != 4){
  print("ERROR: at least start1, end1, start2, end2 are required for input file!")
  q()
}

index_reverse <- which(D_table$start1 > D_table$start2)
D_table[index_reverse, c("start1", "end1", "start2", "end2")] = D_table[index_reverse, c("start2", "end2", "start1", "end1")]

### 色の対象
if(TARGET_column %in% colnames(D_table)){
  scores <- D_table[,TARGET_column]
  if(TARGET_column == "comb"){
    category <- ifelse(scores == "E-P", 1, 2)
  }else{
    bk <- seq(min(scores, na.rm = TRUE), max(scores, na.rm = TRUE), length.out=20)
    category <- cut(D_table[,TARGET_column], breaks=bk, include.lowest = TRUE)
  }
  D_table <- D_table %>% mutate(z=COLOR[as.integer(category)])
}else{
  D_table <- D_table %>% mutate(z=COLOR[1])
}


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

if(FLAG_cairo){
  suppressWarnings(suppressMessages(library(Cairo)))
  CairoPNG(FILE_out, width=10.16, height=2, units="cm", bg="transparent", res=200)
}else{
  png(file=FILE_out, width=10.16, height=2, units="cm", bg="transparent", res=200)
}


par(oma=c(0,0,0,0), mar=c(0,0,0,0))
if(nrow(D_table) > 0){
  Ymax <- D_table %>% mutate(dd=end2-start1) %>% pull(dd) %>% max()/2
}else{
  Ymax = 50
}
if(FLAG_bottom){
  Ymax <- Ymax * -1
}
plot(c(START, END), c(0, Ymax), type='n', xlab="", ylab="", xaxs="i", yaxs="i",bty="n", axes=F)
if(nrow(D_table) > 0){
  dummy <- sapply(1:nrow(D_table), drawEach)
}
dummy <- dev.off()



