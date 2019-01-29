#!/usr/bin/Rscript
# TAD and Compartment summary report

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-d", "--data"), default="NA", help="data directory"),
  make_option(c("--db"), default="NA", help="database file"),
  make_option(c("-o", "--out"), default="NA", help="output directory"),
  make_option(c("-c", "--color"), default="NA", help="colors"),
  make_option(c("--sample"), default="NA", help="sample name separated by ,"),
  make_option(c("--switch"), default="NA", help="sample name of switch pairs, separated by , and : for different combinations ex(A,B:C,D)")
  
)
opt <- parse_args(OptionParser(option_list=option_list))

DIR <- as.character(opt["data"])
SAMPLES <- unlist(strsplit(as.character(opt["sample"]), ","))
COLORS <- unlist(strsplit(as.character(opt["color"]), ","))
SAMPLE_SWITCH <- strsplit(unlist(strsplit(as.character(opt["switch"]), ":")), ",")
DB_sample <- as.character(opt["db"])
DIR_out <- as.character(opt["out"])
if(substring(DIR, nchar(DIR), nchar(DIR)) != "/"){
  DIR <- paste0(DIR, "/")
}
if(substring(DIR_out, nchar(DIR_out), nchar(DIR_out)) != "/"){
  DIR_out <- paste0(DIR_out, "/")
}
DIR_img <- paste0(DIR_out, "img/")
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

options(scipen=10)
options(useFancyQuotes = FALSE)

cut_digits <- function(x, digits=0) {
  return (floor(x * 10^digits) /10^digits)
}

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
  data.frame(Sample=name, length=length, stringsAsFactors = FALSE)
}
D_tadSize <- getTADSize(SAMPLES[1])
for(i in 2:length(SAMPLES)){
  D_tadSize <- rbind(D_tadSize, getTADSize(SAMPLES[i]))
}
p <- ggplot(D_tadSize, aes(x=Sample, y=length/1000)) +geom_jitter(alpha = 0.05, color = "grey40")  +
  geom_boxplot(alpha = 0.3) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))  + theme(legend.position="none") +
  scale_y_log10() +
  annotation_logticks(sides="l") + labs(y="TAD size (kb)", x="")

save_plot(paste0(DIR_out, "img/TAD_size.png"),p, base_height = 5, base_width = 5)

out <- D_tadSize %>% group_by(Sample) %>% summarise(Median=median(length), Average=format(mean(length), digits=2)) %>% as.data.frame()
write.table(out, paste0(DIR_out, "TAD_size.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#=============================================
# Compartment comparison
#=============================================
for(resolution in c("200kb", "40kb")){
  con = dbConnect(SQLite(), DB_sample)
  Query <- paste("select chr, start,", paste0(dQuote(paste0("comp_", SAMPLES)), " as ", dQuote(SAMPLES),  collapse = ", "), " from ", paste0("Comp_", resolution))
  D_comp <- dbGetQuery(con, Query)
  dbDisconnect(con)

  D_comp <- D_comp %>% mutate(id=paste(chr, start, sep=":")) %>% select(-chr, -start)
  D_comp <- D_comp[,c("id", SAMPLES)]
  D_summ <- D_comp %>% tidyr::gather(key = "Sample", value = "compartment", -id)

  out <- D_summ %>% group_by(Sample) %>% summarise(n_A=sum(compartment=="A", na.rm = TRUE), n_B=sum(compartment=="B", na.rm = TRUE)) %>% as.data.frame()
  out2 <- out %>% mutate(p_A=cut_digits(n_A/(n_A+n_B)*100, digits = 2), p_B=cut_digits(n_B/(n_A+n_B)*100, digits = 2)) %>% select(Sample, n_A, p_A, n_B, p_B)
  colnames(out2) <- c("Sample", "n", "%", "n", "%")
  header <- as.data.frame(t(c("", "A", "", "B", "")))
  FILE_out <- paste0(DIR_out, "Compartment_distribution_", resolution, ".txt")
  write.table(header, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  suppressWarnings(write.table(out2, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE))

  D_bar <- out %>% tidyr::gather(key = "compartment", value = "count", -Sample)

  p <- ggplot(D_bar, aes(x = Sample, y = count, fill = compartment)) +scale_fill_manual(values=c("red", "blue"), name="Comp.")+
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
  for(i in 1:length(SAMPLE_SWITCH)){
    con = dbConnect(SQLite(), DB_sample)
    Query <- paste("select chr, start,", paste0(dQuote(paste0("comp_", SAMPLE_SWITCH[[i]])), " as ", dQuote(SAMPLE_SWITCH[[i]]),  collapse = ", "), " from ", paste0("Comp_", resolution))
    D_comp <- dbGetQuery(con, Query)
    dbDisconnect(con)

    num_na <- apply(is.na(D_comp), 1, sum)
    D_comp <- D_comp[num_na==0,]
    D_comp <- D_comp %>% mutate(id=paste(chr, start, sep=":"), comb=paste0(D_comp[,SAMPLE_SWITCH[[i]][1]], D_comp[,SAMPLE_SWITCH[[i]][2]])) %>% select(id, comb)

    out <- D_comp %>% group_by(comb) %>% summarise(n=n()) %>% as.data.frame()
    total <- sum(out$n)
    out2 <- out %>% mutate(p=cut_digits(n/total*100, digits = 1))
    colnames(out2) <- c("comb", "n", "%")
    FILE_out <- paste0(DIR_out, "CompartmentSwitch_", SAMPLE_SWITCH[[i]][1], "_to_", SAMPLE_SWITCH[[i]][2], "_", resolution, ".txt")
    write.table(out2, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat(paste0("Total\t", total, "\n"), file=FILE_out, append = TRUE)

    out <- out %>% mutate(comb=factor(comb, levels = rev(c("AA", "BB", "AB", "BA"))))
    p <- ggplot(out, aes(x = "", y = n, fill = comb)) +scale_fill_manual(values=rev(c("red", "blue", "yellow", "green")), name="Comp.")+
        geom_bar(stat = "identity", position = "fill") +
        scale_y_continuous(labels = percent) +
      labs(x="", y="Compartment (%)")
    save_plot(paste0(DIR_out, "img/CompartmentSwitch_", SAMPLE_SWITCH[[i]][1], "_to_", SAMPLE_SWITCH[[i]][2], "_", resolution, ".png"),p, base_height = 5, base_width = 2.8)
  }
}



#=============================================
# Compartment size comparison
#=============================================
for(resolution in c("200kb", "40kb")){
  con = dbConnect(SQLite(), DB_sample)
  Query <- paste("select chr, start, end, ", paste0(dQuote(paste0("comp_", SAMPLES)), " as ", dQuote(SAMPLES),  collapse = ", "), " from ", paste0("Comp_", resolution))
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

  save_plot(paste0(DIR_out, "img/CompartmentSize_", resolution, ".png"),p, base_height = 5, base_width = 8)


  out.1 <- D_size %>% filter(comp=="A") %>% group_by(sample) %>% summarise(n_A=n(), median_A=median(len), average_A=format(mean(len), digits = 2)) %>% as.data.frame()
  out.2 <- D_size %>% filter(comp=="B") %>% group_by(sample) %>% summarise(n_B=n(), median_B=median(len), average_B=format(mean(len), digits = 2)) %>% as.data.frame()
  out <- dplyr::left_join(out.1, out.2, by="sample")
  colnames(out) <- c("Sample", rep(c("n", "Median", "Average"), 2))

  header <- as.data.frame(t(c("", "Compartment A", "", "", "Compartment B", "", "")))
  FILE_out <- paste0(DIR_out, "CompartmentSize_", resolution, ".txt")
  write.table(header, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  suppressWarnings(write.table(out, FILE_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE))
}



#=============================================
# switch Compartment size comparison
#=============================================
for(resolution in c("200kb", "40kb")){
  for(i in 1:length(SAMPLE_SWITCH)){
    con = dbConnect(SQLite(), DB_sample)
    Query <- paste("select chr, start, end, ", paste0(dQuote(paste0("comp_", SAMPLE_SWITCH[[i]])), " as ", dQuote(SAMPLE_SWITCH[[i]]),  collapse = ", "), " from ", paste0("Comp_", resolution))
    D_comp <- dbGetQuery(con, Query)
    dbDisconnect(con)
  
    num_na <- apply(is.na(D_comp), 1, sum)
    D_comp <- D_comp[num_na==0,]
    D_comp <- D_comp %>% mutate(comb=paste0(D_comp[,SAMPLE_SWITCH[[i]][1]], D_comp[,SAMPLE_SWITCH[[i]][2]])) %>% select(chr, start, end, comb)
    
    
    
    getWidth <- function(ccc){
      DD <- D_comp %>% filter(comb == ccc) %>% select(chr, start, end)
      DD <- makeGRangesFromDataFrame(DD)
      DD <- reduce(DD)
      DD <- data.frame(comb=ccc, len=width(DD), stringsAsFactors = FALSE)
      DD
    }
    D_size <- getWidth("AA")
    for(j in c("BB", "AB", "BA")){
      D_size <- rbind(D_size, getWidth(j))
    }
  
    D_size <- D_size %>% mutate(comb=factor(comb, levels = c("AA", "BB", "AB", "BA")))
    p <- ggplot(D_size, aes(x=comb, y=len/1000)) + geom_jitter(alpha = 0.1, color = "grey40")  +
      geom_boxplot(alpha = 0.3, fill=c("red", "blue", "yellow", "green")) +
      scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),　labels = trans_format("log10", math_format(10^.x)) ) +
      annotation_logticks(sides="l") + labs(y="Compartment size (kb)", x="")
  
    save_plot(paste0(DIR_out, "img/CompartmentSwitchSize_", SAMPLE_SWITCH[[i]][1], "_to_", SAMPLE_SWITCH[[i]][2], "_", resolution, ".png"),p, base_height = 5, base_width = 4)
  
    out <- D_size %>% group_by(comb) %>% summarise(n=n(), Median=median(len), Average=format(mean(len), digits = 2)) %>% as.data.frame()
    write.table(out, paste0(DIR_out, "CompartmentSwitchSize_", SAMPLE_SWITCH[[i]][1], "_to_", SAMPLE_SWITCH[[i]][2], "_", resolution, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}
