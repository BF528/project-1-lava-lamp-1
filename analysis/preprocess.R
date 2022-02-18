#Make sure to load in these packages
BiocManager::install(version = "3.14")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"))
library(sva)
library(affyPLM)
library(affy)
library(AnnotationDbi)
library(hgu133plus2.db)
library(bladderbatch)

###Part 1
#CEL Files were read in and normalized together
Data <- ReadAffy(celfile.path="/projectnb2/bf528/users/lava-lamp/project_1/samples")
Data_LavaLamp_rma <- rma(Data, normalize = TRUE)



#The RLE and NUSE Stats included the median and interquartile range
Data_LavaLamp <- fitPLM(Data, normalize=TRUE, background=TRUE)
RLE(Data_LavaLamp, main="Relative Log Expression Plot")
NUSE(Data_LavaLamp, main="Normalized Unscaled Standard Error Plot")
RLE_Data <- RLE(Data_LavaLamp,type = "stats")
NUSE_Data <- NUSE(Data_LavaLamp,type ="stats")



#The RLE and NUSE Histograms were made. Only the medians were utilized. 
hist(RLE_Data["median",], main="Histogram of the Median RLE", xlab= "Median RLE Values",)
hist(NUSE_Data["median",], main="Histogram of the Median NUSE", xlab= "Median NUSE Values")


###Part 2
#The lavalamp data was then corrected via ComBat
exprs_lavalamp <- exprs(Data_LavaLamp_rma)
lavalamp <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
batch = lavalamp$normalizationcombatbatch
mod = model.matrix(~normalizationcombatmod, data=lavalamp)

#The ComBat package was used 
corrected_lavalamp <- ComBat(dat = exprs_lavalamp, batch = batch, mod = mod)

#The corrected_lavalamp is used to make a CSV file
write.csv(corrected_lavalamp, file = "lavalamp.csv")


###Part 3
#The corrected_lavalamp is transposed
Transpose_corrected_lavalamp <- t(corrected_lavalamp)
Scale_Transpose_corrected_lavalamp <- scale(Transpose_corrected_lavalamp)
Transpose_Again_corrected_lavalamp <- t(Scale_Transpose_corrected_lavalamp)
lavalamp_PCA <- prcomp(Transpose_Again_corrected_lavalamp, scale=FALSE, center=FALSE)
lavalamp_PCA_rotated <- as.data.frame(lavalamp_PCA$rotation)
summary(lavalamp_PCA) #Highest standard deviation observed was 3.9059

library(ggplot2)

#The variance percentages of the data are calculated
lavalamp_PCA_variance <- lavalamp_PCA$sdev^2
lavalamp_PCA_variance_percentage <- lavalamp_PCA_variance/sum(lavalamp_PCA_variance)*100
lavalamp_variance_plot <- barplot (lavalamp_PCA_variance_percentage, 
                          main = "Variance Plot",
                          xlab = "x-axis",
                          ylab = "% Variation")

#The PCA Plot with Outliers was plotted
ggplot(data = lavalamp_PCA_rotated, mapping = aes(x=PC1, y=PC2)) +
  geom_point(color = "purple") +
      labs(title = "The PCA Plot")

#The PC1 and PC2 Boxplots to visualize standard deviation was made
#The greatest outlier was seen around 3 standard deviations away from the mean 
ggplot(data = lavalamp_PCA_rotated, mapping = aes(y=PC2)) +
  geom_boxplot(color = "green") +
      labs(title = "The PC2 BoxPlot")

ggplot(data = lavalamp_PCA_rotated, mapping = aes(y=PC1)) +
  geom_boxplot(color = "red") +
      labs(title = "The PC1 BoxPlot")

#Attempting to remove outliers that are 3 standard deviations away from mean
 which(lavalamp_PCA_rotated$PC1 > mean(lavalamp_PCA_rotated$PC1) + 3*sd(lavalamp_PCA_rotated$PC1)
       |lavalamp_PCA_rotated$PC1 < mean(lavalamp_PCA_rotated$PC1) - 3*sd(lavalamp_PCA_rotated$PC1))
 which(lavalamp_PCA_rotated$PC2 > mean(lavalamp_PCA_rotated$PC2) + 3*sd(lavalamp_PCA_rotated$PC2)
       |lavalamp_PCA_rotated$PC2 < mean(lavalamp_PCA_rotated$PC2) - 3*sd(lavalamp_PCA_rotated$PC2))
 
outlier <- which(!(lavalamp_PCA_rotated$PC1 > mean(lavalamp_PCA_rotated$PC1) + 3*sd(lavalamp_PCA_rotated$PC1)
                   |lavalamp_PCA_rotated$PC1 < mean(lavalamp_PCA_rotated$PC1) - 3*sd(lavalamp_PCA_rotated$PC1)) |
                   lavalamp_PCA_rotated$PC2 > mean(lavalamp_PCA_rotated$PC2) + 3*sd(lavalamp_PCA_rotated$PC2)
                 |lavalamp_PCA_rotated$PC2 < mean(lavalamp_PCA_rotated$PC2) - 3*sd(lavalamp_PCA_rotated$PC2))   

LavaLampFiltered <- Transpose_Again_corrected_lavalamp[,outlier]

lavalampfiltered2 <- prcomp(LavaLampFiltered, scale = FALSE, center = FALSE)

#Final PCA Plot
lavalampfiltered_rotated <- as.data.frame(lavalampfiltered2$rotation)
ggplot(data = lavalampfiltered_rotated, mapping = aes(x=PC1, y=PC2)) +
  geom_point(color = "blue") +
      labs(title = "The Filtered PCA Plot")






