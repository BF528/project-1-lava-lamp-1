library(tidyverse)
#Part 5 Clustering Code 

data <- read.csv(file='triplefiltered.csv', sep = ',')
data$X.1 <- NULL
names(data)[1] <- 'ID'

data_label <- data$ID
data$ID <- NULL

transposed_data <- t(data)

#Scale Data

scaled <- as_tibble(scale(transposed_data))

#Create Distance Matrix

dist_mat <- dist(scaled, method = 'euclidean')

hclust_avg <- hclust(dist_mat, method = 'average')

plot(hclust_avg)

cut_avg <- cutree(hclust_avg, k = 2)



#Heatmap 

hmdata <- read.csv(file='triplefiltered.csv', sep = ',')
hmdata$X.1 <- NULL
names(hmdata)[1] <- 'ID'
genelabel <- hmdata$ID
hmdata$ID <- NULL

hmmatrix <- as.matrix(hmdata)

#Assigning ColSideColors
amatrix <- read.csv(file = '/project/bf528/project_1/doc/proj_metadata.csv')



subtypes <- amatrix$cit.coloncancermolecularsubtype
colcolors <- character(length(subtypes))

for(i in 1:length(subtypes)) {
  if(subtypes[i] == "C3") {
    colcolors[i] = "Red"
  } else {
    colcolors[i] = "Blue"
  }
}

heatmap(hmmatrix, labRow = genelabel, ColSideColors = colcolors)

#Welch T-Test
#Data cleanup, assigning clusters to the expression data and transposing 

ttest <- function(filename, cut_avg) {
  tdata <- read.csv(file=filename, sep = ',')
  tdata$X.1 <- NULL
  names(tdata)[1] <- 'ID'
  tlabel <- tdata$ID
  tdata$ID <- NULL
  transposed <- t(tdata)
  tdf <- data.frame(transposed)
  colnames(tdf) = tlabel
  
  tdf$cluster <- cut_avg
  
  clus1 <- filter(tdf, cluster == 1)
  clus2 <- filter(tdf, cluster == 2)
  clus1$cluster <- NULL
  clus2$cluster <- NULL
  c1label <- colnames(clus1)
  c2label <- colnames(clus2)

  clus1t <- data.frame(clus1)
  clus2t <- data.frame(clus2)
  colnames(clus1t) = c1label
  colnames(clus2t) = c2label
  

  IDs <- character(0)
  tstat <- character(0)
  pval <- character(0)
  
  
  for(i in 1:ncol(clus1t)) {
    t1 <- t.test(clus1t[i], clus2t[i])
    pval[i] = t1$p.value
    tstat[i] = t1$statistic[["t"]]
    IDs[i] = colnames(clus1t[i])
  }
  
  adjp <- p.adjust(pval, method = "fdr")
  
  finaldf <- data.frame(IDs, tstat, pval, adjp)
  
  print("The number of differentially expressed genes at adjp < 0.05 is:")
  print(count(finaldf, adjp<0.05)[2,2])
  sorteddf <- finaldf[order(finaldf$adjp),]  
  return(sorteddf)
}
numsample <- table(cut_avg)

print("Number of samples in the first cluster is:")
print(numsample[["1"]])
print("Number of samples in the second cluster is:")
print(numsample[["2"]])

passalldf <- ttest('triplefiltered.csv', cut_avg)
passchidf <- ttest('chisqfiltered.csv', cut_avg)

write.csv(passalldf, 'welchresultsallfilter.csv')
write.csv(passchidf, 'welchresultchisq.csv')

