# Development of PCPR2 script
# for loops have been removed and code vectorised where possible

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

#sum(apply(X_DataMatrix, 2, var) == 0)

Z_Meta <- meta[samples, c(1:2, 10:13)]

#myPath <- "Documents/PCPR2/" Metabolomics_data <- "X_MetaboMatrix.TXT"
#InterestFactors_data <- "Z_FactorMatrix.TXT"
#Metabo_FilePath = paste(myPath,Metabolomics_data, sep="")
#Factors_FilePath = paste(myPath,InterestFactors_data, sep="")
#X_DataMatrix <- read.delim(Metabo_FilePath, row.names = 1, header = TRUE, sep = "\t")

X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

#Do a PCA to check the distribution of variability explained
brew.pca <- prcomp(X_DataMatrixScaled, scale. = F)
#summary(brew.pca)

# Center the data / Scale the data, edit the parameter "pareto" or "unit" of scaling according to your need
#Z_InterestFactors <- read.delim(Factors_FilePath, sep = "\t", header = TRUE, row.names = 1)
Z_MetaRowN <- nrow(Z_Meta)
Z_MetaColN <- ncol(Z_Meta)
ColNames   <- names(Z_Meta)

# get number of PCs for threshold
pct_threshold <- 0.8 # Amount of variability desired to be explained
sumpca <- summary(brew.pca)$importance
pc_n <- which(sumpca[3, ] >= pct_threshold) %>% min

# Get eigenvectors and eigenvalues
X_DataMatrixScaled_t <- t(X_DataMatrixScaled)
symMat <- X_DataMatrixScaled %*% X_DataMatrixScaled_t
eigenData      <- eigen(symMat)
eigenValues    <- eigenData$values
eigenVecMatrix <- eigenData$vectors

pc_data_matrix <- eigenVecMatrix[, 1:pc_n ]

#Perform linear multiple regression models on each eigenvector with factors of interest as explanatory variables
#Categorical variables should be processed by as.factor, whereas continuous variables should not. 
#To be edited with your factors names

# Convert categorical variables to factors. Put them in varlist
varlist <- c("plate", "caffeine", "bean.type", "roast", "brew.method")
Z_Meta <- Z_Meta %>% mutate_at(vars(varlist), as.factor)

DataCol <- Z_MetaColN +1
# Generate Anova type 3 sums of squares for each component and bind rowwise
# make emtpy matrices
type3mat <- matrix(data = 0, nrow = pc_n, ncol = DataCol ) 
ST_ResidualR2 <- matrix(data = 0, nrow = pc_n, ncol = 2)   

# Run a linear model with the with the eigenvector as the response
type3matrows <- function(f) {
  TotSumSq <- var(f) * (Z_MetaRowN - 1)
  
  # Use all interest factors
  fit <- lm(f ~ ., data = Z_Meta)
  
  # Use type III anova to get the sums of squares for each factor
  # Added singular.ok argument to suppress error message
  AnovaTab   <- Anova(fit, type=3, singular.ok = T)
  SumSq      <- AnovaTab[1] 
  Residuals  <- SumSq[DataCol + 1, ] 
  RR         <- Residuals/TotSumSq
  R2         <- 1 - RR
  ST_ResidualR2Row  <- c(R2, RR)
  type3matRow <- SumSq[1:DataCol + 1, 1]
}

# apply model to n data frames where n = number of PCs
type3mat <- apply(pc_data_matrix, 2, type3matrows) %>% t
colnames(type3mat) <- c(ColNames, "SumSqResiduals")

# Calculate ST_ResidualR2 and give colnames
ST_ResidualR2 <- cbind(1 - (type3mat[, DataCol]), type3mat[, DataCol])
colnames(ST_ResidualR2) <- c("ST_R2", "ST_Residuals")

# Make partial R2 matrix and populate from typeIII matrix
makepartialmat <- function(x) x / (x + type3mat[, DataCol])
partialR2mat   <- apply(type3mat[ , -DataCol], 2, makepartialmat)

# Apply eigenvalues as weights
eigenValues <- eigenData$values[1: pc_n]
weight     <- eigenValues/sum(eigenValues) 
partialR2MatWtProp <- cbind(partialR2mat, ST_ResidualR2[, 1])*weight
colnames(partialR2MatWtProp) <- NULL
pR2Sums <- colSums(partialR2MatWtProp) * 100

bp <- barplot(pR2Sums, xlab = "", ylab = "Weighted Rpartial2", ylim = c(0,70), col = "red", 
              las=2, cex.main = 0.8, main = paste("Original PCPR2 vectorised n =", Z_MetaRowN))
axis(1, at = bp, labels = c(ColNames, "R2"), cex.axis = 0.8, las=2) 
rounded <- round(pR2Sums, 3)
text(bp, pR2Sums, labels = rounded, pos=3, cex = 0.8)
