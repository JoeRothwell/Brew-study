#PCPR2 on either all coffee samples (n=174) or all coffee samples except instants (n=150)
#Final analysis for coffee brew paper

pcpr2 <- function(subset.instants = F) {
  
  library(tidyverse)
  library(MetabolAnalyze)
  library(car)
  
  ints <- read_csv("Brew PT Feb 2015.csv", skip = 8) %>% select(2:201) %>% t
  # "Unknown" was replaced with NA in metadata to avoid aliased coefficients in model
  meta <- read_csv("Brew meta pt order march.csv")
  
  #174 coffees to be included (exclusion of blanks and QCs) or 150 (exclusion of instants)
  samples <- if(subset.instants == F) {
    !is.na(meta$coffee.brew)
  } else {
    !is.na(meta$coffee.brew) & meta$brew.method != "Instant"
  }
  
  #samples <- !is.na(meta$coffee.brew) & meta$brew.method != "Instant"
  
  #features detected in more than 2 samples only
  X_DataMatrix <- ints[samples, ]
  filt <- colSums(X_DataMatrix > 1) > 2
  
  #subset both dimensions. Leaves 3626 and 3572 features for full and subset data.
  X_DataMatrix <- log2(ints[samples, filt])
  print(paste(ncol(X_DataMatrix), "variables in X-matrix"))
  
  #sum(apply(X_DataMatrix, 2, var) == 0)
  Z_Meta <- meta[samples, c(1:2, 10:13)]
  
  # Center the data / Scale the data, edit the parameter "pareto" or "unit"
  X_DataMatrixCentered <- scale(X_DataMatrix, center = TRUE, scale = FALSE)
  X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")
  
  
  
  # get number of rows and colnames
  Z_MetaRowN <- nrow(Z_Meta)
  ColNames <- names(Z_Meta)
  
  # Get eigenvectors and eigenvalues
  X_DataMatrixScaled_t <- t(X_DataMatrixScaled)
  symMat <- X_DataMatrixScaled %*% X_DataMatrixScaled_t
  eigenData      <- eigen(symMat)
  eigenValues    <- eigenData$values
  eigenVecMatrix <- eigenData$vectors
  
  # get number of PCs for threshold
  pct_threshold <- 0.8 # Amount of variability desired to be explained
  pc_n <- which(cumsum(eigenValues)/sum(eigenValues) >= pct_threshold) %>% min
  
  # subset eigenvectors for PCs
  pc_data_matrix <- eigenVecMatrix[, 1:pc_n ]

  # Convert categorical variables to factors. Put them in varlist
  varlist <- c("plate", "caffeine", "bean.type", "roast", "brew.method")
  Z_Meta <- Z_Meta %>% mutate_at(vars(varlist), as.factor)
  
  DataCol <- ncol(Z_Meta) + 1
  
  # Run a linear model with the with the eigenvector as the response
  type3matrows <- function(f) {
    #Edit the linear model with your factors
    #fit <- lm(pc_data_matrix ~ plate + run.order + caffeine + bean.type + 
                #roast + brew.method, data = cbind(pc_data_matrix = f, Z_Meta))
    fit <- lm(pc_data_matrix ~ ., data = Z_Meta)
    # Use type III anova to get the sums of squares for each factor
    # Added singular.ok argument to suppress error message
    AnovaTab   <- Anova(fit, type=3, singular.ok = T)
    SumSq      <- AnovaTab[1] 
    type3matRow <- SumSq[1 : DataCol + 1, 1]
  }
  
  # apply model to each col of eigenvector matrix (to pc_n)
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
  weight     <- eigenValues / sum(eigenValues) 
  partialR2MatWtProp <- cbind(partialR2mat, ST_ResidualR2[, 1]) * weight
  colnames(partialR2MatWtProp) <- NULL
  pR2Sums <- colSums(partialR2MatWtProp) * 100
  
  # Plot data
  bp <- barplot(pR2Sums, ylab = "Weighted Rpartial2", ylim = c(0, 70), xlab = "",
        col = "red", las=2, cex.main = 0.8, main = paste("PCPR2 vectorised mod. n =", Z_MetaRowN))
  axis(1, at = bp, labels = c(ColNames, "R2"), cex.axis = 0.8, las=2) 
  rounded <- round(pR2Sums, 2)
  text(bp, pR2Sums, labels = rounded, pos = 3, cex = 0.8)
  
}

par(mfrow = c(1,2))
pcpr2()
pcpr2(subset.instants = T)
type3mat