# 順序を任意に並び替えたcontact mapを描画する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output png file"),
  make_option(c("--order"),help="order file"),
  make_option(c("--color"), default="matlab", help="color matlab or gentle"),
  make_option(c("--unit"), default="p", help="unit to define score threshold p:percent or v:value"),
  make_option(c("--blur"), default="FALSE", help="if TRUE, make blur image"),
  make_option(c("--min"), default="NULL", help="minimum score for drawing"),
  make_option(c("--max"), default="NULL", help="maximu score for drawing"),
  make_option(c("--restriction"), default="cut", help="cut: only draw bingin, line: draw line"),
  make_option(c("--width"), default="1000", help="width of output figure"),
  make_option(c("--height"), default="NULL", help="height of output figure")
)
opt <- parse_args(OptionParser(option_list=option_list))

c <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                        "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))( 100 )
if(as.character(opt["color"]) == "gentle"){
  suppressWarnings(suppressMessages(library("RColorBrewer")))
  c <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100))  # gentle color
}
if(as.character(opt["color"]) == "red"){
  suppressWarnings(library("RColorBrewer"))
  red <- brewer.pal(9, "Reds")[-1]
  c <- colorRampPalette(red)(100)  # gentle color
}

TakeMiddleV <- function(mat, minVal, maxVal){
  mat.new <- ifelse(mat < minVal, minVal, mat)
  ifelse(mat.new > maxVal, maxVal, mat.new)
}
# %で指定する
TakeMiddleP <- function(mat, minPos, maxPos){
  NUM <- sort(as.numeric(mat), na.last = NA)
  minVal <- NUM[length(NUM)*minPos+1]
  maxVal <- NUM[length(NUM)*maxPos+1]
  TakeMiddleV(mat, minVal, maxVal)
}

Transform <- function(mat){
  d = dim(mat)[1]
  n = dim(mat)[2]
  mat.rev = t(mat[d+1-c(1:d), ])
  mat.rev
}



map <- as.matrix(read.table(as.character(opt["in"])), header=TRUE, check.names = FALSE)
colnames(map) <- rownames(map)
OrderData <- read.table(as.character(opt["order"]), header=FALSE, check.names=FALSE)
rownames(OrderData) <-  OrderData[,1]
# Target <- which(OrderData[,2] > -1)
OrderIndex <- intersect(OrderData[,1], rownames(map))
map <- map[OrderIndex, OrderIndex]

# no centromere and telomere
r <- rownames(map)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Range_cen <- 500000
cen1 <- which(as.character(LocMatrix[,1]) == "I" & (as.numeric(LocMatrix[,2]) >= 3771554 - Range_cen) & (as.numeric(LocMatrix[,3]) <= 3771554 + Range_cen))
cen2 <- which(as.character(LocMatrix[,1]) == "II" & (as.numeric(LocMatrix[,2]) >= 1623505 - Range_cen) & (as.numeric(LocMatrix[,3]) <= 1623505 + Range_cen))
cen3 <- which(as.character(LocMatrix[,1]) == "III" & (as.numeric(LocMatrix[,2]) >= 1103953 - Range_cen) & (as.numeric(LocMatrix[,3]) <= 1103953 + Range_cen))
Range_tel <- 100000
tel_left <- which(as.numeric(LocMatrix[,2]) <= Range_tel)
tel1 <- which(as.character(LocMatrix[,1]) == "I" & (as.numeric(LocMatrix[,3]) >= 5579133 - Range_tel))
tel2 <- which(as.character(LocMatrix[,1]) == "II" & (as.numeric(LocMatrix[,3]) >= 4539804 - Range_tel))
tel3 <- which(as.character(LocMatrix[,1]) == "III" & (as.numeric(LocMatrix[,3]) >= 2452883 - Range_tel))
Exclude_list <- c(cen1, cen2, cen3, tel_left, tel1, tel2, tel3)
map <- map[-Exclude_list, -Exclude_list]


# blur image
if(eval(parse(text=as.character(opt["blur"])))){
  suppressWarnings(suppressMessages(library("spatstat")))
  t <- blur(as.im(map), sigma=.6, bleed=FALSE)
  map <- t$v
}


# binding score > 0の数を求める
Num_binding <- sum(OrderData[rownames(map),2] > 0)
border <- Num_binding / nrow(map)
if(as.character(opt["restriction"]) == "cut"){
  map <- map[1:Num_binding, 1:Num_binding]
}

map <- Transform(map)

if(as.character(opt["unit"]) == "p"){
  if(as.character(opt["min"]) == "NULL"){
    Min <- 0.1
  }else{
    Min <- as.numeric(as.character(opt["min"]))
  }
  if(as.character(opt["max"]) == "NULL"){
    Max <- 0.95
  }else{
    Max <- as.numeric(as.character(opt["max"]))
  }
  map.conv <- TakeMiddleP(map, Min, Max)
}else if(as.character(opt["unit"]) == "v"){
  if(as.character(opt["min"]) == "NULL"){
    Min <- min(map, na.rm=TRUE)
  }else{
    Min <- as.numeric(as.character(opt["min"]))
  }
  if(as.character(opt["max"]) == "NULL"){
    Max <- max(map, na.rm=TRUE)
  }else{
    Max <- as.numeric(as.character(opt["max"]))
  }
  map.conv <- TakeMiddleV(map, Min, Max)
}


c_ex <- c[round(1+(min(map.conv, na.rm=T)-Min)*100/(Max-Min)) : round( (max(map.conv, na.rm=T)-Min)*100/(Max-Min) )]


width <- as.numeric(as.character(opt["width"]))
if(as.character(opt["height"]) == "NULL"){
  height <- width / nrow(map.conv) * ncol(map.conv)
}else{
  width <- as.numeric(as.character(opt["height"]))
}

FILE_OUT <- as.character(opt["out"])
png(file=FILE_OUT, width=width, height=height, units="px", bg="white")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
image(map.conv, col=c, axes=F)
if(as.character(opt["restriction"]) == "line"){
  polygon(c(0,border, border, 0), 1-c(0, 0, border, border), col='black', density=0, lwd=3)
}

dummy <- dev.off()

