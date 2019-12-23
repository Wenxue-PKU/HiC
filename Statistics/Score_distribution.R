# HiC score distribution check

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="matrix file"),
  make_option(c("-o", "--out"),help="output text file"),
  make_option(c("--png"), default="NULL", help="output histogram"),
  make_option(c("-t", "--title"), default="", help="title of graph")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(cowplot)))


FILE_in <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
FILE_png <- as.character(opt["png"])
TITLE <- as.character(opt["title"])

if(grepl(".rds", FILE_in)){
  map <- readRDS(FILE_in)
}else{
  map <- as.matrix(read.table(FILE_in, header=TRUE, check.names = FALSE))
}

### NAは0に置換する
map <- ifelse(is.na(map), 0, map)

D_data <- data.frame(score=as.integer(map))

SCORE_MAX <- max(D_data$score)
CATEGORY <- ifelse(SCORE_MAX > 50, 50, SCORE_MAX)

index_center <- as.integer(nrow(map)/2) + 1
SCORE_center <- map[index_center, index_center]
index_cenRegion <- (index_center-1):(index_center+1)
SCORE_centerR <- mean(map[index_cenRegion, index_cenRegion])
sd <- sd(D_data$score)
med <- median(D_data$score)
ave <- mean(D_data$score)

D_statistics <- data.frame(category=c("center", "cen.area", "sd", "median", "ave."), 
                           value=c(SCORE_center, SCORE_centerR, sd, med, ave))
D_statistics <- D_statistics %>% mutate(value=format(value, digits = 3))



D_table <- rbind(D_data %>% filter(score != 0) %>% mutate(cate_sameNum=ntile(score, CATEGORY-1) + 1),
                 D_data %>% filter(score == 0) %>% mutate(cate_sameNum=1))
D_table <- D_table %>% mutate(cate_sameWidth=cut_interval(score, n = CATEGORY))
GRAPH_data <- D_statistics %>% mutate(out=paste(category, value, sep=":")) %>% pull(out)
GRAPH_data <- paste0(GRAPH_data[1], ",  ", GRAPH_data[2], "\n", paste(GRAPH_data[3:5], collapse = ", "))

D_summay <- rbind(
  D_table %>% group_by(cate_sameNum) %>% summarize(min=min(score), max=max(score), n=n(), p=format(n()/nrow(D_table) * 100, digits = 3)) %>%
    as.data.frame() %>% mutate(out=paste0(min, "-", max, ": ", n, " (", p, "%)")) %>% select(out),
  
  D_table %>% group_by(cate_sameWidth) %>% summarize(min=min(score), max=max(score), n=n(), p=format(n()/nrow(D_table) * 100, digits = 3)) %>%
    as.data.frame() %>% mutate(out=paste0(min, "-", max, ": ", n, " (", p, "%)")) %>% select(out)
)

write.table(D_summay, FILE_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


if(FILE_png != "NULL"){
  D_graph <- D_table %>% mutate(score=ifelse(score > CATEGORY, CATEGORY, score)) %>% select(score)
  D_graph <- D_graph %>% group_by(score) %>% summarize(n=n()) %>% as.data.frame()
    
    
  p <- ggplot(D_graph, aes(x=score, y=n)) +  geom_bar(stat = "identity") +
    theme_bw() + theme(text = element_text(size=12), legend.position="top") +
    labs(x="", y="", title=TITLE, subtitle = GRAPH_data)
  
  save_plot(FILE_png, p, base_width = 6, base_height = 4)
  
}


