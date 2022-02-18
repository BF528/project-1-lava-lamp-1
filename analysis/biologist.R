if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.14")
BiocManager::install("hgu133plus2.db")

if (!require("biomaRt", quietly = TRUE)){
  install.packages("biomaRt", quietly = TRUE)
}
BiocManager::install("biomaRt", force = TRUE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")

library(GSEABase)
library(GO.db)

library(tidyverse)
library(BiocManager)
library(biomaRt) #needed for ensembl
library(readr)
library(dplyr)

load_expression <- function(filepath) {
  read_data <- read_delim(filepath, delim = ',')
  colnames(read_data) <- c('probeids', 't', 'p', 'padj')
  #separate(data = read_data, col = 'p', into = c('p', 'padj'), sep = ',')
  #tib_data <- as_tibble(read_data)
  return (read_data)
}

expr <- load_expression('data/differential_expression_results.csv')



require(hgu133plus2.db)
columns(hgu133plus2.db)
keys = pull(expr['probeids'])
#keys

match <- select(hgu133plus2.db, keys, columns = c('PROBEID', 'SYMBOL'))
match <- match[!(is.na(match$SYMBOL) | match$SYMBOL==""), ]
match$SYMBOL[duplicated(match$SYMBOL)] <- ''
match <- match[!(match$SYMBOL==''),]

reduce_data <- function(df){

}

expr <- expr %>%
  filter(probeids %in% match$PROBEID)
expr$Symbol <- match$SYMBOL

match <- match %>%
  filter(PROBEID %in% expr$probeids)

#select the top 1000 up- and down-regulated
chisq <- read_delim('data/welchresultchisq.csv', delim = ',')


columns(hgu133plus2.db)
chisq_symb <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(chisq$IDs), 
                                    columns = ("SYMBOL"))

chisq_dist_symb <- chisq_symb[!duplicated(chisq_symb$PROBEID),]

chisq_filt <- merge(chisq_dist_symb, chisq, by.x = 'PROBEID', by.y = 'IDs' )

#remove duplicates and keep one probe per sample based on adjusted pvalue
chisq_filt <- chisq_filt %>%
  group_by(SYMBOL) %>% filter(adjp == min(adjp)) %>% ungroup(SYMBOL)

chisq_filt <- chisq_filt[-c(1,3)]
  


#Top 1000 up and down regulated
up_reg <- top_n(chisq_filt, 1000, tstat)
down_reg <- top_n(chisq_filt, -1000, tstat)

#top 10 up and down regulated
top_ten <- top_n(up_reg, 10, pval)
bot_ten <- top_n(down_reg, -10, pval)

write.csv(top_ten, 'Top10upregGenes.csv')
write.csv(bot_ten, 'Top10downregGenes.csv')

#fisher test

#contingency table

cont_table <- function(gs, ngs, none){  #gs = in gene set, ngs = not in gs(top1000list)
  exp_in_gs <- length(intersect(gs,ngs))
  nexp_in_gs <- length(intersect(ngs, none))
  exp_in_nogs <- length(ngs) - exp_in_gs
  nexp_in_nogs <- length(none) - nexp_in_gs
  return(c(exp_in_gs, exp_in_nogs, nexp_in_gs, nexp_in_nogs))
}

up_no_exp <- subset(chisq_filt, !chisq_filt$SYMBOL %in% up_reg$SYMBOL)
down_no_exp <- subset(chisq_filt, !chisq_filt$SYMBOL %in% down_reg$SYMBOL)

#The gene sets used were from the Hallmark, KEGG, and GO databases

#Hallmark

hallmark <- getGmt('data/h.all.v7.5.1.symbols.gmt')
length(names(hallmark)) #number of gene sets considered


#result dataframes
hall_stats <- data.frame(gs_name = character(), stat_est = as.numeric(), pval = as.numeric())

#hallmark fisher test
for (i in 1:length(hallmark)){
  hall_id <- geneIds(hallmark[i])
  
  hall_up_test <- fisher.test(matrix(cont_table(up_reg $ SYMBOL, 
                                      hall_id[[names(hall_id)]], up_no_exp$SYMBOL),nrow=2))
  hall_down_test <- fisher.test(matrix(cont_table(down_reg $ SYMBOL, 
                                        hall_id[[names(hall_id)]], down_no_exp$SYMBOL), nrow=2))
  hall_stats[nrow(hall_stats) + 1,] <- c(names(hall_id), hall_up_test$estimate, hall_up_test$p.value)
  hall_stats[nrow(hall_stats) + 1,] <- c(names(hall_id), hall_down_test$estimate, hall_down_test$p.value)
  
}

hall_stats
hall_stats$adjp <- p.adjust(hall_stats$pval, method = "BH",
                          n = length(hall_stats$pval))


hall_up_test
hall_down_test

sig_hall_stats = hall_stats %>% filter(pval <0.05)
count(sig_hall_stats) #number of significantly enriched genes
hall_3 <- top_n(hall_stats, 3, stat_est)
write.csv(hall_3, 'halltop3.csv')


#Kegg

kegg <- getGmt('data/c2.cp.kegg.v7.5.1.symbols.gmt')
length(names(kegg))

kegg_stats <- data.frame(gs_name = character(), stat_est = as.numeric(), pval = as.numeric())

#Kegg fisher test
for (i in 1:length(kegg))
{
  keggid <- geneIds(kegg[i])
  
  kegg_up_test <- fisher.test(matrix(cont_table(up_reg $ SYMBOL, 
                                                keggid[[names(keggid)]], up_no_exp$SYMBOL),nrow=2))
  kegg_down_test <- fisher.test(matrix(cont_table(down_reg $ SYMBOL, 
                                                  keggid[[names(keggid)]], down_no_exp$SYMBOL), nrow=2))
  kegg_stats[nrow(kegg_stats) + 1,] <- c(names(keggid), kegg_up_test$estimate, kegg_up_test$p.value)
  kegg_stats[nrow(kegg_stats) + 1,] <- c(names(keggid), kegg_down_test$estimate, kegg_down_test$p.value)
  
  }
kegg_stats
kegg_stats$adjp <- p.adjust(kegg_stats$pval, method = "BH",
                            n = length(kegg_stats$pval))

kegg_up_test
kegg_down_test

sig_kegg_stats = kegg_stats %>% filter(pval <0.05)

count(sig_kegg_stats)
kegg3 <- top_n(kegg_stats, 3, stat_est)
write.csv(kegg3, 'keggtop3.csv')


#GO

go <- getGmt('data/c5.go.v7.5.1.symbols.gmt')
length(names(go))

go_stats <- data.frame(gs_name = character(), stat_est = as.numeric(), pval = as.numeric())


for (i in 1:length(go))
{
  go_id <- geneIds(go[i])
  
  go_up_test <- fisher.test(matrix(cont_table(up_reg $ SYMBOL, 
                                                go_id[[names(go_id)]], up_no_exp$SYMBOL),nrow=2))
  go_down_test <- fisher.test(matrix(cont_table(down_reg $ SYMBOL, go_id[[names(go_id)]], down_no_exp$SYMBOL), nrow=2))

  go_stats[nrow(go_stats) + 1,] <- c(names(go_id), go_up_test$estimate, go_up_test$p.value)
  go_stats[nrow(go_stats) + 1,] <- c(names(go_id), go_down_test$estimate, go_down_test$p.value)
}

go_stats
go_stats$adjp <- p.adjust(go_stats$pval, method = "BH",
                          n = length(go_stats$pval))

go_up_test
go_down_test

sig_go_stats = go_stats %>% filter(pval < 0.05)

count(sig_go_stats)
go3 <- top_n(go_stats, 3, stat_est)
write.csv(go3, 'gotop3.csv')

write.csv(sig_hall_stats,'sighall.csv')
write.csv(sig_kegg_stats,'sigkegg.csv')
write.csv(sig_go_stats,'siggo.csv')
