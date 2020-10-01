#!/usr/bin/Rscript
# Distance curve from matrix for intra compartment

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix or rds file(s)"),
  make_option(c("-o", "--out"),default="NA", help="output text file"),
  make_option(c("--normalize"), default="Probability", help="Average will be 1, Probability: score were divided by total read, NA:without normalization"),
  make_option(c("--compartment"), default="NA", help="compartment definition, provide file with chr, start, end, compartment, domain_id")

)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(pbapply)))
suppressWarnings(suppressMessages(library(GenomicRanges)))

Normalization <- as.character(opt["normalize"])
FILE_comp <- as.character(opt["compartment"])
FILE_out <- as.character(opt["out"])
FILE_in <- unlist(strsplit(as.character(opt["in"]), ","))

if(FILE_comp != "NA"){
  D_comp <- fread(FILE_comp, header=T) %>% as.data.frame()
  G_comp <- GRanges(D_comp)
}

getAverageDistance <- function(file){
  if(file.exists(file)){
    map <- readRDS(file)
  }else{
    file <- sub(".matrix", ".rds", file)
    if(file.exists(file)){
      map <- as.matrix(fread(file, header=TRUE, check.names = FALSE))
    }else{
      cat(file, " not exists\n")
      q()
    }
  }
  
  ### Normalization
  if(Normalization == "Average"){
    d <- nrow(map)
    map <- map / sum(map, na.rm=TRUE) * d * d
  }else if(Normalization == "Probability"){
    map <- map / sum(map, na.rm=TRUE)
  }else if(Normalization == "NA"){
    # don't change anything
  }else{
    cat("Normalization parameter is wrong")
    q()
  }
  
  ### compartment information
  D_Location <- data.frame(id=rownames(map)) %>% mutate(str=id) %>% tidyr::separate(str, c("chr", "start", "end"), ":")
  G_Location <- GRanges(D_Location)
  ov <- findOverlaps(G_Location, G_comp)
  overlaps_len <- width(pintersect(G_Location[queryHits(ov)], G_comp[subjectHits(ov)]))
  D_Location <- data.frame(id=G_Location$id[queryHits(ov)], compartment=G_comp$compartment[subjectHits(ov)], 
                            domain_id=G_comp$domain_id[subjectHits(ov)], overlaps_len, start=start(G_Location)[queryHits(ov)])
  D_Location <- D_Location %>% group_by(id) %>% filter(overlaps_len==max(overlaps_len)) %>% as.data.frame()
  rm(G_Location)
  
  ### domain_idごとにスコアを抽出
  getScore <- function(i){
    df <- D_Location %>% filter(domain_id==i)
    id_list <- df %>% pull(id)
    if(length(id_list) > 2){
      id_combination <- combn(id_list,2)
      df2 <- data.frame(id1=id_combination[1,], id2=id_combination[2,])
      df2 <- dplyr::left_join(df2, df %>% select(id, start) %>% dplyr::rename(loc1=start), by=c("id1"="id"))
      df2 <- dplyr::left_join(df2, df %>% select(id, start) %>% dplyr::rename(loc2=start), by=c("id2"="id"))
      df2 <- df2 %>% mutate(dis=loc2-loc1)
      scores <- map[df2 %>% select(id1, id2) %>% as.matrix()]
      df <- data.frame(scores, distance = df2 %>% pull(dis))
      df2 <- df %>% filter(!is.na(scores)) %>% group_by(distance) %>% summarize(n=n(), sum=sum(scores)) %>% as.data.frame()
      df2 <- df2 %>% mutate(compartment=D_Location %>% filter(domain_id==i) %>% head(n=1) %>% pull(compartment))
      df2
    }else{
      return(NULL)
    }
  }
  D_count <- do.call(rbind, pblapply(D_Location %>% filter(compartment != "N") %>% pull(domain_id), getScore))
  D_count <- D_count %>% group_by(distance, compartment) %>% summarize(n=sum(n), sum=sum(sum)) %>% as.data.frame()
  D_count
}

D_table <- do.call(rbind, pblapply(FILE_in, getAverageDistance))
D_table <- D_table %>% group_by(distance, compartment) %>% summarize(n=sum(n), sum=sum(sum)) %>% as.data.frame()
D_table <- D_table %>% mutate(average=sum/n) %>% select(-n, -sum)
write.table(D_table, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


