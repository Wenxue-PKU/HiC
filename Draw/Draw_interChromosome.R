# inter-chromosomeのマップ

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("--organism"), default="human", help="organism"),
  make_option(c("--no_chrY"), default="TRUE", help="remove chrY"),
  make_option(c("--enzyme"), default="length", help="normalize by chromosome length or enzyme name (MboI or HindIII)"),
  make_option(c("--max"), default="1.5", help="max log2 score"),
  make_option(c("--min"), default="-1.5", help="min log2 score"),
  make_option(c("--name"), default="", help="title of graph"),
  make_option(c("-o", "--out"),help="output png file")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_matrix <- "X:/hideki_projects/2017-08-08_HiC_human_senescent/data/CR_Sen/InterChromosome.matrix"
ORGANISM <- "human"
MAX <- 1.5
MIN <- -1.5
NAME <- "CR_Sen"

FILE_matrix <- as.character(opt["in"])
ORGANISM <- as.character(opt["organism"])
MAX <- as.numeric(as.character(opt["max"]))
MIN <- as.numeric(as.character(opt["min"]))
NAME <- as.character(opt["name"])
ENZYME <- as.character(opt["enzyme"])

suppressWarnings(suppressMessages(library("RColorBrewer")))
c <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100))  # gentle color


if(ORGANISM == "human"){
  chr_names <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  if(ENZYME == "length"){
    chr_length <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  }else if(ENZYME == "HindIII"){
    chr_length <- c(64394, 71482, 59635, 58924, 55088, 51310, 45516, 43186, 34642, 38015, 39325, 38986, 29411, 25751, 22932, 19902, 18938, 22770, 11541, 15704, 10088, 7581, 45073, 7429)
  }else if(ENZYME == "MboI"){
    chr_length <- c(576356, 578289, 471036, 439279, 427116, 406041, 389012, 347577, 299461, 331194, 325829, 331816, 223049, 219601, 210286, 214856, 219287, 178582, 171356, 155186, 83815, 96365, 372192, 60004)
  }
}


  

if(eval(parse(text=as.character(opt["no_chrY"])))){
  index_chrY <- which(chr_names == "chrY")
  chr_names <- chr_names[-index_chrY]
  chr_length <- chr_length[-index_chrY]
}


map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
map <- map[chr_names, chr_names]
map[cbind(chr_names, chr_names)] <- NA

expected <- matrix(chr_length, nrow=length(chr_names), ncol=length(chr_names), byrow = TRUE) * 
  matrix(chr_length, nrow=length(chr_names), ncol=length(chr_names), byrow = FALSE)
colnames(expected) <- colnames(map)
rownames(expected) <- rownames(map)


map <- map / sum(map, na.rm=TRUE) * length(map) * length(map)
expected <- expected/ sum(expected, na.rm=TRUE) * length(expected) * length(expected)
obs_exp <- log2(map / expected)

Transform <- function(mat){
  d = dim(mat)[1]
  n = dim(mat)[2]
  mat.rev = t(mat[d+1-c(1:d), ])
  mat.rev
}

obs_exp <- ifelse(obs_exp > MAX, MAX, obs_exp)
obs_exp <- ifelse(obs_exp < MIN, MIN, obs_exp)

c_ex <- c[round(1+(min(obs_exp, na.rm=T)-MIN)*100/(MAX-MIN)) : round( (max(obs_exp, na.rm=T)-MIN)*100/(MAX-MIN) )]

FILE_OUT <- as.character(opt["out"])
png(file=FILE_OUT, width=500, height=500, units="px", bg="white")
par(oma=c(0,0,2,0), mar=c(1,4,4,1))
image(Transform(obs_exp), col=c_ex, axes=F)
axis(side=3,at=seq(0,1, length=length(chr_names)),labels=chr_names, las=2)
axis(side=2,at=seq(0,1, length=length(chr_names)),labels=rev(chr_names), las=2)
mtext(NAME, outer = TRUE, side=3, cex=1.2, line=0)
