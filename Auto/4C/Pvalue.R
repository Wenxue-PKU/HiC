#!/usr/bin/Rscript
# P-value calculation

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="bed graph file"),
  make_option(c("-o", "--out"),help="p-value list")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))

DATA <- read.table(as.character(opt["in"]), header=F, sep="\t", stringsAsFactors = FALSE, check.names = FALSE, col.names = c("chr", "start","end", "score"))

# EBV以外の場所のaverage
average <- DATA %>% filter(chr != "EBV") %>% pull(score) %>% mean()

p <- ppois(DATA$score, lambda = average, lower.tail = F)

logp <- log10(p)*-1
logp <- ifelse(is.infinite(logp), 1000, logp)

DATA <- DATA %>% mutate(logp=logp) %>% arrange(chr, start) %>% select(chr, start, end, logp)

write.table(DATA, as.character(opt["out"]), quote = F, sep="\t", eol = "\n", row.names = F, col.names = F)

