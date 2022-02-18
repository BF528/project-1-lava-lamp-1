library(tidyverse)
path = 'preprocess.csv'

data <- read.csv(file = path)

##Filter out probes where 20% samples are < log2(15)
data$F1 <- FALSE

for(i in 2:nrow(data)) {
  count = 0 
  for(k in 2:(ncol(data)-1)) {
    if((data[i,k]) > log2(15)) {
      count <- count + 1 
    }
  }
  if(count > ncol(data) * 0.2) {
    data[i, ncol(data)] <- TRUE
  }
}

f1data <- filter(data, data['F1'] == TRUE)
f1data$F1 <- NULL



##Variance significantly Different from all probe sets 
f1data$var <- apply(f1data[,-1],1,var)
totalmed <- median(f1data$var)

qst <- qchisq(0.01,133)
#f1data <- cbind(f1data, qst )


teststatistic <- function(x) { 
  N = 34
  s = totalmed
  t = (N-1)*(s)/x
  return(t)
}

data_var <- f1data[, ncol(f1data), drop=FALSE]
data_var$tst <- apply(data_var, 1, teststatistic)
f1data$tst <- data_var$tst

f2data <- filter(f1data, f1data['tst'] < qst)




#Calculating cov and filtering for the third time 
f2data$cov <- apply(f2data[2:(ncol(f2data)-2)],1,function(x) sd(x)/mean(x))
filtered <- filter(f2data, f2data['cov'] > 0.186)

print("The number which passed all the first filters was: ") 
print(count(f1data))
print("The number which passed all the first filters was: ") 
print(count(f2data))
print("The number which passed all three filters was: ") 
print(count(filtered))

filtered$var <- NULL
filtered$tst <- NULL
filtered$cov <- NULL
f2data$var <- NULL
f2data$tst <- NULL
f2data$cov <- NULL
#End of part 4 
write.csv(filtered, 'triplefiltered.csv')
write.csv(f2data, 'chisqfiltered.csv')
