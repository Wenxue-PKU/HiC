# making mask for inter-chromosome, intra-chromosomal inter-chromosomal

RESOLUTION <- "10kb"
FILE_object <- paste("Z:/Data/2017-03-28_HiC_pombe_revising_paper/WT-36C2/", RESOLUTION, "/ICN/ALL.rds", sep="")
FILE_interChr <- paste("Z:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/mask/Mask_interChr_", RESOLUTION, ".rds", sep="")
FILE_interArm <- paste("Z:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/mask/Mask_interArm_", RESOLUTION, ".rds", sep="")
FILE_intraChr <- paste("Z:/Project/2017-03-28_HiC_pombe_revising_paper/out/2017-05-25_all_statistics/mask/Mask_intraChr_", RESOLUTION, ".rds", sep="")

map_ref <- readRDS(FILE_object)

r <- rownames(map_ref)
LocList <- strsplit(r, ":")
LocMatrix <- matrix(unlist(LocList), ncol=3, byrow=TRUE)
Resolution <- as.numeric(LocMatrix[1,3]) - as.numeric(LocMatrix[1,2]) + 1
NUM_LINE <- nrow(map_ref)

#=======================================
# inter-chromosome
#=======================================
# interChr <- function(x){
#   as.character(LocMatrix[x,1]) != as.character(LocMatrix[,1])
# }
# mask.interChr <- cbind(sapply(1:NUM_LINE, interChr))
# rownames(mask.interChr) <- r
# colnames(mask.interChr) <- r
# 
# saveRDS(mask.interChr, FILE_interChr)


#=======================================
# intra-chromosome
#=======================================
intraChr <- function(x){
  as.character(LocMatrix[x,1]) == as.character(LocMatrix[,1])
}
mask.intraChr <- cbind(sapply(1:NUM_LINE, intraChr))
rownames(mask.intraChr) <- r
colnames(mask.intraChr) <- r

saveRDS(mask.intraChr, FILE_intraChr)





#=======================================
# intra-chromosome inter-arm
#=======================================
# ChromSection <- rep(0, nrow(LocMatrix))
# CHR1 <- as.character(LocMatrix[,1]) =="I"
# CHR2 <- as.character(LocMatrix[,1]) =="II"
# CHR3 <- as.character(LocMatrix[,1]) =="III"
# ChromSection[which(CHR1 & as.numeric(LocMatrix[,2]) < 3753687)] <- 10
# ChromSection[which(CHR1 & as.numeric(LocMatrix[,3]) > 3789421)] <- 11
# ChromSection[which(CHR2 & as.numeric(LocMatrix[,2]) < 1602264)] <- 20
# ChromSection[which(CHR2 & as.numeric(LocMatrix[,3]) > 1644747)] <- 21
# ChromSection[which(CHR3 & as.numeric(LocMatrix[,2]) < 1070904)] <- 30
# ChromSection[which(CHR3 & as.numeric(LocMatrix[,3]) > 1137003)] <- 31
# interArm <- function(x){
#   abs(ChromSection[x] - ChromSection) == 1
# }
# mask.interArm <- cbind(sapply(1:NUM_LINE, interArm))
# rownames(mask.interArm) <- r
# colnames(mask.interArm) <- r
# 
# saveRDS(mask.interArm, FILE_interArm)


