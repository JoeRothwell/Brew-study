# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
library(MetabolAnalyze)
library(car)
library(tidyverse)
# load(file="brewdata_log2_pl.Rda")

# Explaining variables (Y) -----------------------------------------------------------------------------------------
meta <- read_csv("Brew meta pt order march.csv")

#get the vector for all or non-instant coffees only (removing QC and blank)
all.cof  <- !is.na(meta$coffee.brew)
non.inst <- meta$brew.method != "Instant" & !is.na(meta$coffee.brew)


#subset the metadata for coffee samples only and variables of interest only
Y_InterestFactors <- meta[all.cof, c(1:2, 10:13)]
#non-instant coffees only
#Y_InterestFactors <- meta[non.inst, c(1:2, 10:13)]
dim(Y_InterestFactors)

Y_InterestFactorsRowN <- nrow(Y_InterestFactors)
Y_InterestFactorsColN <- ncol(Y_InterestFactors)
ColNames <- names(Y_InterestFactors)

# Variables to explain (X) (untargeted metabolomics data) ----------------------------------------------------------

ints <- read_csv("Brew PT Feb 2015.csv", skip=8) %>% select(2:201) %>% t
ints[ints == 0] <- 1
# Subset for selected observations (coffee samples) only
ints <- log2(ints[all.cof, ])
#ints <- log2(ints[non.inst, ])
dim(ints)

#test for zero variance columns
filt <- which(apply(ints, 2, var) == 0)
#one only, 3667

#metabolomics data for all or non-instant coffees only
X_DataMatrix <- ints[ , -filt]
#test for NAs and get dimensions
anyNA(X_DataMatrix)
dim(X_DataMatrix)
#dimensions are 174, 3669

#X_DataMatrix <- outcome
#X_DataMatrix$sample_name = NULL

# Center the data / Scale the data, edit the parameter "pareto" or "unit" of scaling according to your need
X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "unit")

#antibug
stopifnot(!anyNA(X_DataMatrixScaled))

time_scale <- proc.time() # Initialization of the time
#X_DataMatrixScaled = scale(X_DataMatrix, center=T, scale=T)
cat("total scale :", proc.time()[3]-time_scale[3], "s \n")

# Obtain eigenvectors

pct_threshold = .8 # Amount of variability desired to be explained, to be edited with your preferences


time_all <- proc.time() # Initialization of the time

# PCPR2 function -------------------------------------------------------------------------
X_DataMatrixScaled_transposed = t(X_DataMatrixScaled)
Mat2      <- X_DataMatrixScaled %*% X_DataMatrixScaled_transposed
eigenData <- eigen(Mat2)
eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues / eigenValuesSum
my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1) {
  my_sum_2  = my_sum_2 - percents_PCs[i]
  if ((my_sum_2) <= pct_threshold) {
    my_counter_2 = my_counter_2 + 1
  }
}
if (my_counter_2 < 3) {
  pc_n  = 3
} else {
  pc_n = my_counter_2
}
pc_data_matrix <-
  matrix(
    data = 0,
    nrow = (Y_InterestFactorsRowN * pc_n),
    ncol = 1
  )
#it = 0
for (i in 1:pc_n) {
  for (j in 1:Y_InterestFactorsRowN) {
    it <- (i - 1) * Y_InterestFactorsRowN + j
    pc_data_matrix[it, 1] = eigenVectorsMatrix[j, i]
  }
}
AAA <- Y_InterestFactors[rep(1:Y_InterestFactorsRowN, pc_n),]
names(AAA) <- colnames(Y_InterestFactors)
Data <- cbind(pc_data_matrix, AAA)

# Perform linear multiple regression models on each eigenvector with factors of interest as explanatory variables
#Categorical variables should be processed by as.factor, whereas continuous variables should not. To be edited with your factors names
names(Data)

#for (i in 1:dim(Data)[2]) {
#  if(nlevel(Data[,i])>1){
#    Data[,i] = as.factor(Data[,i])
#  }
#}

Data$plate <- as.factor(Data$plate)
#Data$sample.type <- as.factor(Data$sample.type)
#Data$coffee.brew <- as.factor(Data$coffee.brew)
Data$caffeine <- as.factor(Data$caffeine)
Data$bean.type <- as.factor(Data$bean.type)
Data$roast <- as.factor(Data$roast)
Data$brew.method <- as.factor(Data$brew.method)

#check all are factors
str(Data)


DataCol <- ncol(Data)
typeIIIMatrix <- matrix(data = 0, nrow = (pc_n), ncol = (DataCol) )
ST_ResidualR2 <- matrix(data = 0, nrow = (pc_n), ncol = 2)

for (i in 1:pc_n){
  y = (((i-1) * Y_InterestFactorsRowN) + 1)
  TotSumSq <- var(Data[y:(((i-1) * Y_InterestFactorsRowN) + Y_InterestFactorsRowN),1]) * (Y_InterestFactorsRowN-1)
  #Edit the linear model with your factors
  data = Data[y:(((i-1)*Y_InterestFactorsRowN) + Y_InterestFactorsRowN),]
  f <- as.formula(paste(names(data)[1], "~", paste(names(data)[-1], collapse=" + ")))
  Model <- lm((f), data) # data[,2:dim(data)[2]]
  AnalysisVariance <- Anova(Model, type=c(3))
  SumSq <- AnalysisVariance[1] 
  Residuals <- SumSq[DataCol + 1, ]
  RR <- Residuals/TotSumSq 
  R2 = 1 - RR
  ST_ResidualR2[i,] <- c(R2,RR)
  ST_ResidualR2_Names <- c("ST_R2", "ST_Residuals")
  colnames(ST_ResidualR2) = ST_ResidualR2_Names
  for (j in 1:(DataCol)){
    typeIIIMatrix[i,j] = as.numeric(SumSq[j+1,1])
    typeIIIMatrixNames <- c(ColNames, "SumSqResiduals")
    colnames(typeIIIMatrix) = typeIIIMatrixNames
  }
}
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
pR2Sums <-colSums(partialR2MatrixWtProp)*100
plotnames = c( ColNames, "R2")
# x11()
bp <- barplot(pR2Sums,  xlab = "Factors", ylab = "Weighted Rpartial2", ylim=c(0,70),col = c("red"), las=2)
axis(1, at = bp, labels = plotnames, xlab = "Factors", cex.axis = 0.45, las=1)
values = pR2Sums
new_values = round(values, 3)
text(bp, pR2Sums, labels = new_values, pos=3, cex = 0.8)

cat("total time :", proc.time()[3]-time_all[3], "s \n")

# bp <- barplot(pR2Sums,  xlab = "Factors", ylab = "Weighted Rpartial2", ylim= c(0,60),col = c("red"), las=2)
# axis(1, at = bp, labels = plotnames, xlab = "Factors", cex.axis = 0.45, las=1)
# values = pR2Sums
# new_values = round(values, 3)
# text(bp,pR2Sums,labels = new_values, pos=3, cex = 0.8)

x11()
test = pR2Sums
names(test) = c( ColNames, "R2")
test = sort(test)
bp <- barplot(test,  xlab = "", ylab = "Weighted Rpartial2", ylim= c(0,50), col = c("red"),
              las=2, cex.names = 0.70)
axis(1, at = bp,labels =rep("",length(test)),  xlab = "Factors", cex.axis = 0.45, las=1)
text(bp,test,labels = round(test,1), pos=3, cex = 0.85)
