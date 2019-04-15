#!/usr/bin/Rscript
# Create intra-chromosome sliding window matrices from database

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--db"), default="NA", help="fragment db"),
  make_option(c("-r", "--resolution"), default="NA", help="resolution (bp)"),
  make_option(c("-s", "--sliding"), default="NA", help="sliding window (bp)"),
  make_option(c("-o", "--out"), default="NA", help="matrices name xxx.matrix"),
  make_option(c("-c", "--chr"), default="NA", help="target chromosome")
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(dplyr)))
options(useFancyQuotes = FALSE)

#=============================================
# test data
#=============================================
# DB_frag <- "X:/hideki_projects/416_20181226_HiC_tempera_EBV/data/EBV_LCL_olap_fragment.db"
# SLIDING <- 1000
# RESOLUTION <- 5000
# TARGET_CHR <- "EBV"


DB_frag <- as.character(opt["db"])
RESOLUTION <- as.numeric(as.character(opt["resolution"]))
SLIDING <- as.numeric(as.character(opt["sliding"]))
FILE_matrix <- as.character(opt["out"])
FILE_object <- sub(".matrix", ".rds", FILE_matrix)
TARGET_CHR <- as.character(opt["chr"])
SURROUNDING_BIN <- ((RESOLUTION / SLIDING) - 1)/2



con = dbConnect(SQLite(), DB_frag)
query <- paste0("select start1, end1, start2, end2, score from fragment where chr1==chr2 and chr2==", dQuote(TARGET_CHR), " and fragNum1 != fragNum2;")
D_table <- dbGetQuery(con, query)
dbDisconnect(con)
rm(con)

conv2bin <- function(s, e){
  as.integer((s+e)/2/SLIDING)*SLIDING
}

D_table <- D_table %>% mutate(distance=abs(start1+end1-start2-end2)/2)
D_table <- D_table %>% mutate(score=ifelse(distance < 10000, score*2, score))
D_table <- D_table %>% select(-distance) 
D_table <- D_table %>% mutate(bin1=conv2bin(D_table$start1, D_table$end1), bin2=conv2bin(D_table$start2, D_table$end2))
D_table <- D_table %>% select(-start1, -end1, -start2, -end2)
D_table <- D_table %>% group_by(bin1, bin2) %>% summarize(score=sum(score)) %>% as.data.frame()
BIN_max <- D_table %>% pull(bin2) %>% max()

BINs <- seq(0, BIN_max, by=SLIDING)

D_table <- D_table %>% mutate(bin1=factor(bin1, levels = BINs, ordered = TRUE), bin2=factor(bin2, levels = BINs, ordered = TRUE))
map <- D_table %>% tidyr::spread(key=bin2,value=score, fill=0, drop=FALSE)
# この時点では、mapは右上部分しかスコアがないことに注意すること
rm(D_table)
map <- map[,-1]
map <- as.matrix(map)
colnames(map) <- BINs
rownames(map) <- BINs

### 下半分にコピー
for(i in 1:(nrow(map)-1)){
  map[cbind((i+1):nrow(map), i)] <- map[cbind(i, (i+1):nrow(map))]
}

if(SURROUNDING_BIN != 0){
  map_new <- map
  for (i in (1:nrow(map))){
    CalcAverage <- function(m){
      sum(map[max(1,i-SURROUNDING_BIN):min(nrow(map),i+SURROUNDING_BIN),
               max(1,m-SURROUNDING_BIN):min(ncol(map),m+SURROUNDING_BIN)], na.rm = TRUE)
    }
    
    score <- sapply(i:nrow(map), CalcAverage)
    map_new[cbind(i, i:nrow(map))] <- score
    map_new[cbind(i:nrow(map), i)] <- score
  }
  map <- map_new
  rm(map_new)
}


BIN_names <- paste(TARGET_CHR, BINs, BINs + RESOLUTION - 1, sep=":")
colnames(map) <- BIN_names
rownames(map) <- BIN_names

saveRDS(map, FILE_object)
write.table(map, FILE_matrix, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)








