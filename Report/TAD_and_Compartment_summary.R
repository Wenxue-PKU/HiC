#!/usr/bin/Rscript
# TAD and Compartment summary report

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--data"), default="NA", help="data directory"),
  make_option(c("--db"), default="NA", help="database file"),
  make_option(c("-o", "--out"), default="NA", help="output directory"),
  make_option(c("-c", "--color"), default="NA", help="colors"),
  make_option(c("--sample"), default="NA", help="sample name separated by ,"),
  make_option(c("--switch"), default="NA", help="sample name of switch pairs, separated by ,")
  
)
opt <- parse_args(OptionParser(option_list=option_list))

DIR <- as.character(opt["data"])
SAMPLES <- unlist(strsplit(as.character(opt["sample"]), ","))
COLORS <- unlist(strsplit(as.character(opt["color"]), ","))
SAMPLE_SWITCH <- unlist(strsplit(as.character(opt["switch"]), ","))
DB_sample <- as.character(opt["db"])
DIR_out <- as.character(opt["out"])
DIR_img <- paste0(DIR_out, "img/")

if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
  DIR <- paste0(DIR, "/")
}
if(substring(DIR_out, nchar(DIR_out), nchar(DIR_out)) != "/"){
  DIR_out <- paste0(DIR_out, "/")
}

if(!dir.exists(DIR_img)){
  dir.create(DIR_img)
}


suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(GenomicRanges)))

#=============================================
# Register TAD border strength to database
#=============================================
getTAD <- function(name) {
  file <- paste0(DIR, name, "/TAD_40kb_score.txt")
  DATA_head <- read.table(file, header=FALSE, nrows = 5, stringsAsFactors = FALSE)
  classes <- sapply(DATA_head, class)
  classes[-c(1,2,3,5,6)] <- "NULL"
  DD <- read.table(file, header=FALSE, colClasses = classes, stringsAsFactors = FALSE)
  rm(DATA_head)
  colnames(DD) <- c("chr", "start", "end", paste0("tadSt_", name), paste0("tad_", name))
  DD
}
D_tad <- getTAD(SAMPLES[1])
for(i in 2:length(SAMPLES)){
  D_tad <- dplyr::left_join(D_tad, getTAD(SAMPLES[i]), by=c("chr", "start", "end"))
}
con = dbConnect(SQLite(), DB_sample)
dbWriteTable(con, "TAD", D_tad, row.names= FALSE, overwrite=TRUE)
dbDisconnect(con)

#=============================================
# Register Compartment to database
#=============================================
getCompartment <- function(name, resolution){
  file <- paste0(DIR, name, "/Compartment_", resolution, ".txt")
  DATA_head <- read.table(file, header=FALSE, nrows = 500, stringsAsFactors = FALSE)
  classes <- sapply(DATA_head, class)
  if(length(classes) == 6){
    classes[5] <- "NULL"
  }
  DD <- read.table(file, header=FALSE, colClasses = classes, stringsAsFactors = FALSE)
  rm(DATA_head)
  colnames(DD) <- c("chr", "start", "end", paste0("pca_", name), paste0("comp_", name))
  DD
}
for(rr in c("200kb", "40kb")){
  D_comp <- getCompartment(SAMPLES[1], rr)
  for(i in 2:length(SAMPLES)){
    D_comp <- dplyr::left_join(D_comp, getCompartment(SAMPLES[i], rr), by=c("chr", "start", "end"))
  }
  con = dbConnect(SQLite(), DB_sample)
  dbWriteTable(con, paste0("Comp_", rr), D_comp, row.names= FALSE, overwrite=TRUE)
  dbDisconnect(con)
}

#=============================================
# TAD size comparison
#=============================================
getTADSize <- function(name) {
  file <- paste0(DIR, name, "/TAD_40kb.txt")
  DATA_head <- read.table(file, header=FALSE, nrows = 5, stringsAsFactors = FALSE)
  classes <- sapply(DATA_head, class)
  DD <- read.table(file, header=FALSE, colClasses = classes, stringsAsFactors = FALSE)
  rm(DATA_head)
  colnames(DD) <- c("chr", "start", "end")
  length <- DD[,"end"] - DD[,"start"]
  data.frame(sample=name, length=length, stringsAsFactors = FALSE)
}
D_tadSize <- getTADSize(SAMPLES[1])
for(i in 2:length(SAMPLES)){
  D_tadSize <- rbind(D_tadSize, getTADSize(SAMPLES[i]))
}
p <- ggplot(D_tadSize, aes(x=sample, y=length/1000)) +geom_jitter(alpha = 0.05, color = "grey40")  +
  geom_boxplot(alpha = 0.3) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position="none") +
  scale_y_log10() +
  annotation_logticks(sides="l") + labs(y="TAD size (kb)", x="")

save_plot(paste0(DIR_out, "img/TAD_size.png"),p, base_height = 5, base_width = 5)

out <- D_tadSize %>% group_by(sample) %>% summarise(median=median(length), mean=mean(length)) %>% as.data.frame()
write.table(out, paste0(DIR_out, "TAD_size.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#=============================================
# Compartment comparison
#=============================================
for(resolution in c("200kb", "40kb")){
  con = dbConnect(SQLite(), DB_sample)
  Query <- paste("select chr, start,", paste0("comp_", SAMPLES, " as ", SAMPLES,  collapse = ", "), " from ", paste0("Comp_", resolution))
  D_comp <- dbGetQuery(con, Query)
  dbDisconnect(con)

  D_comp <- D_comp %>% mutate(id=paste(chr, start, sep=":")) %>% select(-chr, -start)
  D_comp <- D_comp[,c("id", SAMPLES)]
  D_summ <- D_comp %>% tidyr::gather(key = "sample", value = "compartment", -id)

  out <- D_summ %>% group_by(sample) %>% summarise(A=sum(compartment=="A", na.rm = TRUE), B=sum(compartment=="B", na.rm = TRUE)) %>% as.data.frame()
  write.table(out, paste0(DIR_out, "compartment_distribution_", resolution, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  D_bar <- out %>% tidyr::gather(key = "compartment", value = "count", -sample)

  p <- ggplot(D_bar, aes(x = sample, y = count, fill = compartment)) +scale_fill_manual(values=c("red", "blue"), name="Comp.")+
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = percent) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x="", y="Compartment (%)")
  save_plot(paste0(DIR_out, "img/Compartment_AB_ratio", resolution, ".png"),p, base_height = 5, base_width = 7)
}

#=============================================
# Compartment change
#=============================================
for(resolution in c("200kb", "40kb")){
  con = dbConnect(SQLite(), DB_sample)
  Query <- paste("select chr, start,", paste0("comp_", SAMPLE_SWITCH, " as ", SAMPLE_SWITCH,  collapse = ", "), " from ", paste0("Comp_", resolution))
  D_comp <- dbGetQuery(con, Query)
  dbDisconnect(con)
  
  num_na <- apply(is.na(D_comp), 1, sum)
  D_comp <- D_comp[num_na==0,]
  D_comp <- D_comp %>% mutate(id=paste(chr, start, sep=":"), comb=paste0(D_comp[,SAMPLE_SWITCH[1]], D_comp[,SAMPLE_SWITCH[2]])) %>% select(id, comb)

  out <- D_comp %>% group_by(comb) %>% summarise(count=n()) %>% as.data.frame()
  write.table(out, paste0(DIR_out, "compartment_switch_", resolution, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  out <- out %>% mutate(comb=factor(comb, levels = rev(c("AA", "BB", "AB", "BA"))))
  p <- ggplot(out, aes(x = "", y = count, fill = comb)) +scale_fill_manual(values=rev(c("red", "blue", "yellow", "green")), name="Comp.")+
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = percent) +
    labs(x="", y="Compartment (%)")
  save_plot(paste0(DIR_out, "img/Compartment_switch_", resolution, ".png"),p, base_height = 5, base_width = 2.8)
}


#=============================================
# Compartment size comparison
#=============================================
for(resolution in c("200kb", "40kb")){
  con = dbConnect(SQLite(), DB_sample)
  Query <- paste("select chr, start, end, ", paste0("comp_", SAMPLES, " as ", SAMPLES,  collapse = ", "), " from ", paste0("Comp_", resolution))
  D_comp <- dbGetQuery(con, Query)
  dbDisconnect(con)
  
  
  getComp_size <- function(name){
    A <- D_comp %>% filter(!!as.name(name) == "A") %>% select(chr, start, end)
    A <- makeGRangesFromDataFrame(A)
    A <- reduce(A)
    A <- data.frame(sample=name, comp="A", len=width(A), stringsAsFactors = FALSE)
    
    B <- D_comp %>% filter(!!as.name(name) == "B") %>% select(chr, start, end)
    B <- makeGRangesFromDataFrame(B)
    B <- reduce(B)
    B <- data.frame(sample=name, comp="B", len=width(B), stringsAsFactors = FALSE)
    
    rbind(A, B)
  }
  D_size <- getComp_size(SAMPLES[1])
  for(i in 2:length(SAMPLES)){
    D_size <- rbind(D_size, getComp_size(SAMPLES[i]))
  }
  
  p <- ggplot(D_size, aes(x=comp, y=len/1000)) + facet_grid(. ~ sample) + geom_jitter(alpha = 0.1, color = "grey40")  +
    geom_boxplot(alpha = 0.3, fill=rep(c("red", "blue"), times=length(SAMPLES))) +
    scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),　labels = trans_format("log10", math_format(10^.x)) ) +
    annotation_logticks(sides="l") + labs(y="Compartment size (kb)", x="")

  save_plot(paste0(DIR_out, "img/Compartment_size_", resolution, ".png"),p, base_height = 5, base_width = 8)
  
  
  out.1 <- D_size %>% filter(comp=="A") %>% group_by(sample) %>% summarise(n_A=n(), median_A=median(len), average_A=mean(len)) %>% as.data.frame()
  out.2 <- D_size %>% filter(comp=="B") %>% group_by(sample) %>% summarise(n_B=n(), median_B=median(len), average_B=mean(len)) %>% as.data.frame()
  
  out <- dplyr::left_join(out.1, out.2, by="sample")
  write.table(out, paste0(DIR_out, "Compartment_size_", resolution, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}

#=============================================
# switch Compartment size comparison
#=============================================
for(resolution in c("200kb", "40kb")){
  con = dbConnect(SQLite(), DB_sample)
  Query <- paste("select chr, start, end, ", paste0("comp_", SAMPLE_SWITCH, " as ", SAMPLE_SWITCH,  collapse = ", "), " from ", paste0("Comp_", resolution))
  D_comp <- dbGetQuery(con, Query)
  dbDisconnect(con)

  num_na <- apply(is.na(D_comp), 1, sum)
  D_comp <- D_comp[num_na==0,]
  D_comp <- D_comp %>% mutate(comb=paste0(D_comp[,SAMPLE_SWITCH[1]], D_comp[,SAMPLE_SWITCH[2]])) %>% select(chr, start, end, comb)
  
  
  
  getWidth <- function(ccc){
    DD <- D_comp %>% filter(comb == ccc) %>% select(chr, start, end)
    DD <- makeGRangesFromDataFrame(DD)
    DD <- reduce(DD)
    DD <- data.frame(comb=ccc, len=width(DD), stringsAsFactors = FALSE)
    DD
  }
  D_size <- getWidth("AA")
  for(i in c("BB", "AB", "BA")){
    D_size <- rbind(D_size, getWidth(i))
  }

  D_size <- D_size %>% mutate(comb=factor(comb, levels = c("AA", "BB", "AB", "BA")))
  p <- ggplot(D_size, aes(x=comb, y=len/1000)) + geom_jitter(alpha = 0.1, color = "grey40")  +
    geom_boxplot(alpha = 0.3, fill=c("red", "blue", "yellow", "green")) +
    scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),　labels = trans_format("log10", math_format(10^.x)) ) +
    annotation_logticks(sides="l") + labs(y="Compartment size (kb)", x="")

  save_plot(paste0(DIR_out, "img/Compartment_switch_size_", resolution, ".png"),p, base_height = 5, base_width = 4)

  out <- D_size %>% group_by(comb) %>% summarise(n=n(), median=median(len), average=mean(len)) %>% as.data.frame()
  write.table(out, paste0(DIR_out, "Compartment_switch_size_", resolution, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
