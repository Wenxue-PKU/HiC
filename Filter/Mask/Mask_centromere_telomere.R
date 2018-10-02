# making mask for centromere and telomere

RESOLUTION <- "10kb"
FILE_object <- paste("Z:/Data/2017-03-28_HiC_pombe_revising_paper/WT-36C2/", RESOLUTION, "/ICN/ALL.rds", sep="")
FILE_telomere <- paste("Z:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/mask/Mask_telomere_", RESOLUTION, ".rds", sep="")
FILE_centromere <- paste("Z:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/mask/Mask_centromere_", RESOLUTION, ".rds", sep="")

map_ref <- readRDS(FILE_object)

r <- rownames(map_ref)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
NUM_LINE <- nrow(map_ref)



#=======================================
# telomere
#=======================================
ChromSection <- rep(0, nrow(LocMatrix))
LENGTH_tel <- 500000

CHR1 <- as.character(LocMatrix[,1]) =="I"
CHR2 <- as.character(LocMatrix[,1]) =="II"
CHR3 <- as.character(LocMatrix[,1]) =="III"
ChromSection[which(CHR1 & as.numeric(LocMatrix[,2]) < LENGTH_tel)] <- 1
ChromSection[which(CHR2 & as.numeric(LocMatrix[,2]) < LENGTH_tel)] <- 1
ChromSection[which(CHR1 & as.numeric(LocMatrix[,3]) > 5579133 - LENGTH_tel)] <- 1
ChromSection[which(CHR2 & as.numeric(LocMatrix[,3]) > 4539804 - LENGTH_tel)] <- 1

telomereAso <- function(x){
  ChromSection[x]==1 & ChromSection == 1
}
mask.tel <- cbind(sapply(1:NUM_LINE, telomereAso))
rownames(mask.tel) <- r
colnames(mask.tel) <- r

saveRDS(mask.tel, FILE_telomere )



#=======================================
# centromere
#=======================================
ChromSection <- rep(0, nrow(LocMatrix))
LENGTH_cen <- 500000
cen1 <- (3753687 + 3789421)/2
cen2 <- (1602264 + 1644747)/2
cen3 <- (1070904 + 1137003)/2
CHR1 <- as.character(LocMatrix[,1]) =="I"
CHR2 <- as.character(LocMatrix[,1]) =="II"
CHR3 <- as.character(LocMatrix[,1]) =="III"
ChromSection[which(CHR1 & as.numeric(LocMatrix[,3]) > cen1 - LENGTH_cen & as.numeric(LocMatrix[,2]) < cen1 + LENGTH_cen)] <- 1
ChromSection[which(CHR2 & as.numeric(LocMatrix[,3]) > cen2 - LENGTH_cen & as.numeric(LocMatrix[,2]) < cen2 + LENGTH_cen)] <- 1
ChromSection[which(CHR3 & as.numeric(LocMatrix[,3]) > cen3 - LENGTH_cen & as.numeric(LocMatrix[,2]) < cen3 + LENGTH_cen)] <- 1

cenAso <- function(x){
  ChromSection[x]==1 & ChromSection == 1 & as.character(LocMatrix[x,1]) != as.character(LocMatrix[,1])
}
mask.cen <- cbind(sapply(1:NUM_LINE, cenAso))
rownames(mask.cen) <- r
colnames(mask.cen) <- r

saveRDS(mask.cen, FILE_centromere )



