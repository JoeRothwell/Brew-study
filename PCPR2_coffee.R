# PCPR2 on either all coffee samples (n=174) or all coffee samples except instants (n=150)
# Final analysis for coffee brew paper

# Original PCPR2 ----

# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
library(MetabolAnalyze)
library(car)
library(tidyverse)

ints <- read_csv("Brew PT Feb 2015.csv", skip = 8) %>% select(2:201) %>% t
#"Unknown" was replaced with NA in metadata to avoid aliased coefficients in model
meta <- read_csv("Brew meta pt order march.csv")

#174 coffees to be included (exclusion of blanks and QCs) or 150 (exclusion of instants)
#samples <- !is.na(meta$coffee.brew)
samples <- !is.na(meta$coffee.brew) & meta$brew.method != "Instant"

#features detected in more than 2 samples only
X_DataMatrix <- ints[samples, ]
filt <- colSums(X_DataMatrix > 1) > 2

#subset both dimensions. Leaves 3626 and 3572 features for full and subset data.
X_DataMatrix <- log2(ints[samples, filt])
#X_DataMatrix <- log2(ints)[samples, -3667]

#sum(apply(X_DataMatrix, 2, var) == 0)

Z_InterestFactors <- meta[samples, c(1:2, 10:13)]

#myPath <- "Documents/PCPR2/" Metabolomics_data <- "X_MetaboMatrix.TXT"
#InterestFactors_data <- "Z_FactorMatrix.TXT"
#Metabo_FilePath = paste(myPath,Metabolomics_data, sep="")
#Factors_FilePath = paste(myPath,InterestFactors_data, sep="")
#X_DataMatrix <- read.delim(Metabo_FilePath, row.names = 1, header = TRUE, sep = "\t")

X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

#Do a PCA to check the distribution of variability explained
#brew.pca <- prcomp(X_DataMatrixScaled, scale. = F)
#summary(brew.pca)

# Center the data / Scale the data, edit the parameter "pareto" or "unit" of scaling according to your need
#Z_InterestFactors <- read.delim(Factors_FilePath, sep = "\t", header = TRUE, row.names = 1)
Z_InterestFactorsRowN <- nrow(Z_InterestFactors)
Z_InterestFactorsColN <- ncol(Z_InterestFactors)
ColNames <- names(Z_InterestFactors)

# Obtain eigenvectors
pct_threshold = .8 # Amount of variability desired to be explained, to be edited with your preferences
X_DataMatrixScaled_transposed = t(X_DataMatrixScaled)
Mat2 <- X_DataMatrixScaled %*% X_DataMatrixScaled_transposed
eigenData <- eigen(Mat2)
eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues /eigenValuesSum
my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1){
  my_sum_2 = my_sum_2 - percents_PCs[i]
  if ((my_sum_2) <= pct_threshold ){
    my_counter_2 = my_counter_2 + 1
  }
}
if (my_counter_2 < 3){ pc_n =3
}else {
  pc_n = my_counter_2
}
pc_data_matrix <- matrix(data = 0, nrow = (Z_InterestFactorsRowN*pc_n), ncol = 1)
mycounter = 0
for (i in 1:pc_n){
  for (j in 1:Z_InterestFactorsRowN){
    mycounter <- mycounter + 1                    
    pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
  }
}
AAA <- Z_InterestFactors[rep(1:Z_InterestFactorsRowN,pc_n),]
Data <- cbind(pc_data_matrix, AAA)

#Perform linear multiple regression models on each eigenvector with factors of interest as explanatory variables
#Categorical variables should be processed by as.factor, whereas continuous variables should not. 
#To be edited with your factors names

Z_InterestFactors$plate <- as.factor(Z_InterestFactors$plate)
#Z_InterestFactors$sample.type <- as.factor(Z_InterestFactors$sample.type)
#Z_InterestFactors$coffee.brew <- as.factor(Z_InterestFactors$coffee.brew)
Z_InterestFactors$caffeine <- as.factor(Z_InterestFactors$caffeine)
Z_InterestFactors$bean.type <- as.factor(Z_InterestFactors$bean.type)
Z_InterestFactors$roast<- as.factor(Z_InterestFactors$roast)
Z_InterestFactors$brew.method <- as.factor(Z_InterestFactors$brew.method)

DataCol <- ncol(Data)
typeIIIMatrix <- matrix(data = 0, nrow = (pc_n), ncol = (DataCol) ) 
ST_ResidualR2 <- matrix(data = 0, nrow = (pc_n), ncol = 2)  

# Run type III ANOVA on each PC
for (i in 1:pc_n){                                        
  y = (((i-1)*Z_InterestFactorsRowN)+1)                                                                                                                      
  TotSumSq <- var(Data[y:(((i- 1)*Z_InterestFactorsRowN)+Z_InterestFactorsRowN),1])*(Z_InterestFactorsRowN-1)
  #Edit the linear model with your factors
  Model <- lm((pc_data_matrix ~ plate + run.order + caffeine + bean.type + 
                 roast + 
                 brew.method), 
              Data[y:(((i-1) * Z_InterestFactorsRowN) + Z_InterestFactorsRowN), ])
  #Note: added singular.ok argument to suppress error message
  AnalysisVariance <- Anova(Model, type=c(3), singular.ok = T)
  SumSq     <- AnalysisVariance[1] 
  Residuals <- SumSq[DataCol + 1, ] 
  RR <- Residuals/TotSumSq
  R2 = 1 - RR
  ST_ResidualR2[i,]   <- c(R2, RR)
  ST_ResidualR2_Names <- c("ST_R2", "ST_Residuals")
  colnames(ST_ResidualR2) = ST_ResidualR2_Names
  for (j in 1:(DataCol)){
    typeIIIMatrix[i,j] = as.numeric(SumSq[j + 1, 1])
    typeIIIMatrixNames <- c(ColNames, "SumSqResiduals")
    colnames(typeIIIMatrix) = typeIIIMatrixNames
  } 
}

#Create partial R2 matrix
partialR2Matrix <- matrix(data = 0, nrow = (pc_n), ncol = (DataCol-1) )
for (i in 1:pc_n){
  for (j in 1:(DataCol-1)){
    partialR2Matrix[i,j] = typeIIIMatrix[i,j] / (typeIIIMatrix[i,(DataCol)] + typeIIIMatrix[i,j]) 
  }
}

partialR2MatrixWtProp <- matrix(data = 0, nrow = (pc_n), ncol = (DataCol)) 
for (i in 1:pc_n){
  weight = eigenValues[i]/sum(eigenValues[1:pc_n]) 
  for (j in 1:DataCol-1){
    partialR2MatrixWtProp[i,j] = partialR2Matrix[i,j]*weight
    partialR2MatrixWtProp[i,DataCol] = ST_ResidualR2[i,1]*weight 
  }
}

pR2Sums <- colSums(partialR2MatrixWtProp)*100 
plotnames = c( ColNames, "R2")
bp <- barplot(pR2Sums, xlab = "", ylab = "Weighted Rpartial2", ylim = c(0,70), col = c("red"), 
              las=2, main = paste("Original PCPR2 n =", ev_n))
axis(1, at = bp, labels = plotnames, cex.axis = 0.8, las=2) 
values = pR2Sums
new_values = round(values, 3)
text(bp, pR2Sums, labels = new_values, pos=3, cex = 0.8)

output <- data.frame(plotnames, pR2Sums)

# Simplified PCPR2 ----
# (code tidied up a bit)

#PCPR2 on either all coffee samples (n=174) or all coffee samples except instants (n=150)
#Final analysis for coffee brew paper

# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
library(MetabolAnalyze)
library(car)
library(tidyverse)

ints <- read_csv("Brew PT Feb 2015.csv", skip = 8) %>% select(2:201) %>% t
#"Unknown" was replaced with NA in metadata to avoid aliased coefficients in model
meta <- read_csv("Brew meta pt order march.csv")

#174 coffees to be included (exclusion of blanks and QCs) or 150 (exclusion of instants)
#samples <- !is.na(meta$coffee.brew)
samples <- !is.na(meta$coffee.brew) & meta$brew.method != "Instant"

#features detected in more than 2 samples only
X_DataMatrix <- ints[samples, ]
filt <- colSums(X_DataMatrix > 1) > 2

#subset both dimensions. Leaves 3626 and 3572 features for full and subset data.
X_DataMatrix <- log2(ints[samples, filt])
#X_DataMatrix <- log2(ints)[samples, -3667]

#sum(apply(X_DataMatrix, 2, var) == 0)

Z_Meta <- meta[samples, c(1:2, 10:13)]

X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

#saveRDS(X_DataMatrixScaled, "Brew_data_PCPR2.rds")
X_DataMatrixScaled <- readRDS("Brew_data_PCPR2.rds")

#Do a PCA to check the distribution of variability explained
#brew.pca <- prcomp(X_DataMatrixScaled, scale. = F)
#summary(brew.pca)

# Center the data / Scale the data, edit the parameter "pareto" or "unit"
Z_MetaRowN <- nrow(Z_Meta)
Z_MetaColN <- ncol(Z_Meta)
ColNames   <- names(Z_Meta)

# Obtain eigenvectors
pct_threshold = .8 # Amount of variability desired to be explained
X_DataMatrixScaled_transposed = t(X_DataMatrixScaled)
Mat2 <- X_DataMatrixScaled %*% X_DataMatrixScaled_transposed
eigenData   <- eigen(Mat2)
eigenValues <- eigenData$values
ev_n               <- length(eigenValues)
eigenVectorsMatrix <- eigenData$vectors

percents_PCs       <- eigenValues / sum(eigenValues)
my_counter_2 <- 0
my_sum_2 <- 1

# get number of PCs for threshold
for (i in ev_n:1){
  my_sum_2 = my_sum_2 - percents_PCs[i]
  if (my_sum_2 <= pct_threshold ) my_counter_2 <- my_counter_2 + 1
}

if (my_counter_2 < 3) pc_n <- 3 else pc_n <- my_counter_2

# stack the eigenvectors in one column
pc_data_matrix <- matrix(data = 0, nrow = Z_MetaRowN * pc_n, ncol = 1)
mycounter <- 0
for (i in 1:pc_n){
  for (j in 1:Z_MetaRowN){
    mycounter <- mycounter + 1                    
    pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, i]
  }
}

# Repeat metadata vertically n times for PCs and bind to eigenvectors
AAA <- Z_Meta[rep(1 : Z_MetaRowN,pc_n), ]
Data <- cbind(pc_data_matrix, AAA)

#Perform linear multiple regression models on each eigenvector with factors of interest as explanatory variables
#Categorical variables should be processed by as.factor, whereas continuous variables should not. 
#To be edited with your factors names

Z_Meta$plate <- as.factor(Z_Meta$plate)
#Z_Meta$sample.type <- as.factor(Z_Meta$sample.type)
#Z_Meta$coffee.brew <- as.factor(Z_Meta$coffee.brew)
Z_Meta$caffeine <- as.factor(Z_Meta$caffeine)
Z_Meta$bean.type <- as.factor(Z_Meta$bean.type)
Z_Meta$roast <- as.factor(Z_Meta$roast)
Z_Meta$brew.method <- as.factor(Z_Meta$brew.method)


DataCol <- ncol(Data)
# Generate Anova type 3 sums of squares for each component and bind rowwise
# make emtpy matrices
type3mat <- matrix(data = 0, nrow = pc_n, ncol = DataCol ) 
ST_ResidualR2 <- matrix(data = 0, nrow = pc_n, ncol = 2)  

# Run type III ANOVA on each PC (i)
for (i in 1 : pc_n){ 
  
  y <- ((i-1) * Z_MetaRowN) + 1
  yy <- (i-1) * Z_MetaRowN
  
  # Total sums of squares 
  TotSumSq <- var(Data[y : (yy + Z_MetaRowN), 1]) * (Z_MetaRowN - 1)
  #Edit the linear model with your factors
  Model <- lm(pc_data_matrix ~ plate + run.order + caffeine + bean.type + 
                roast + brew.method, Data[y:(yy + Z_MetaRowN), ])
  
  AnalysisVariance <- Anova(Model, type=3, singular.ok = T) #added singular.ok to suppress error
  SumSq            <- AnalysisVariance[1] 
  
  # Populate typeIII matrix from sums of squares for each factor and give colnames
  for (j in 1:DataCol) {
    type3mat[i, j]     <- as.numeric(SumSq[j + 1, 1])
    colnames(type3mat) <- c(ColNames, "SumSqResiduals")
  } 
  
  # Populate ST_ResidualR2 and give colnames
  Residuals <- SumSq[DataCol + 1, ] 
  RR <- Residuals/TotSumSq
  R2 <- 1 - RR
  ST_ResidualR2[i, ]      <- c(R2, RR)
  colnames(ST_ResidualR2) <- c("ST_R2", "ST_Residuals")
}

#Create partial R2 matrix and populate from typeIII matrix
partialR2Matrix <- matrix(data = 0, nrow = pc_n, ncol = DataCol-1 )
for (i in 1:pc_n){
  for (j in 1 : DataCol-1){
    partialR2Matrix[i,j] = type3mat[i,j] / (type3mat[i, DataCol] + type3mat[i,j]) 
  }
}

# Created weighted matrix from previous
partialR2MatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = DataCol) 
for (i in 1:pc_n){
  weight <- eigenValues[i]/sum(eigenValues[1 : pc_n]) 
  for (j in 1 : DataCol-1){
    partialR2MatrixWtProp[i, j] <- partialR2Matrix[i,j]*weight
    partialR2MatrixWtProp[i, DataCol] <- ST_ResidualR2[i,1]*weight 
  }
}

pR2Sums <- colSums(partialR2MatrixWtProp)*100 

bp <- barplot(pR2Sums, xlab = "", ylab = "Weighted Rpartial2", ylim = c(0,70), 
              col = "red", las=2, main = paste("Original PCPR2 simp. n =", ev_n))
axis(1, at = bp, labels = c(ColNames, "R2"), cex.axis = 0.8, las=2) 
rounded <- round(pR2Sums, 3)
text(bp, pR2Sums, labels = rounded, pos = 3, cex = 0.8)


output <- data.frame(Factor = c(ColNames, "R2"), pR2Sums)
