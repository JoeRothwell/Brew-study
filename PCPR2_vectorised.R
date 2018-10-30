#PCPR2 on either all coffee samples (n=174) or all coffee samples except instants (n=150)
#Final analysis for coffee brew paper

pcpr2 <- function(subset.instants = F) {
  
  library(tidyverse)
  library(MetabolAnalyze)
  library(car)
  
  ints <- read_csv("Brew PT Feb 2015.csv", skip = 8) %>% select(2:201) %>% t
  #"Unknown" was replaced with NA in metadata to avoid aliased coefficients in model
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
  
  #sum(apply(X_DataMatrix, 2, var) == 0)
  Z_Meta <- meta[samples, c(1:2, 10:13)]
  
  # Center the data / Scale the data, edit the parameter "pareto" or "unit"
  X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
  X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")
  
  # Do a PCA to check the distribution of variability explained
  brew.pca <- prcomp(X_DataMatrixScaled, retx = F, scale. = F)
  
  Z_MetaRowN <- nrow(Z_Meta)
  Z_MetaColN <- ncol(Z_Meta)
  ColNames <- names(Z_Meta)
  
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
  
  # Replicate interest factors data and bind rowwise
  #AAA <- Z_Meta[rep(1 : Z_MetaRowN, pc_n), ]
  
  # Stack eigenvectors and bind to interest factors
  #pc_data_matrix <- c(eigenVecMatrix[, 1:pc_n ])
  pc_data_matrix <- eigenVecMatrix[, 1:pc_n ]
  #Data <- cbind(pc_data_matrix, AAA)
  
  # Make a dataframe for each PC and put in list
  #dflist <- apply(eigenVecMatrix[, 1:pc_n], 2, function(x) cbind(pc_data_matrix = x, Z_Meta))
  
  # Convert categorical variables to factors. Put them in varlist
  varlist <- c("plate", "caffeine", "bean.type", "roast", "brew.method")
  Z_Meta <- Z_Meta %>% mutate_at(vars(varlist), as.factor)
  
  #DataCol <- ncol(dflist[[1]])
  DataCol <- Z_MetaColN +1
  
  # Run a linear model with the with the eigenvector as the response
  type3matrows <- function(f) {
    TotSumSq <- var(f) * (Z_MetaRowN - 1)
    
    #Edit the linear model with your factors
    fit <- lm(f ~ plate + run.order + caffeine + bean.type + 
                  roast + brew.method, 
              data = Z_Meta)
    
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
  #ll <- lapply(dflist, type3matrows)
  #type3mat           <- do.call(rbind, ll) 
  colnames(type3mat) <- c(ColNames, "SumSqResiduals")
  #type3mat[, DataCol]
  
  
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

  # Plot data
  bp <- barplot(pR2Sums, ylab = "Weighted Rpartial2", ylim = c(0, 70), xlab = "",
      col = "red", las=2, cex.main = 0.8, main = paste("PCPR2 vectorised NEW n =", Z_MetaRowN))
  axis(1, at = bp, labels = c(ColNames, "R2"), cex.axis = 0.8, las=2) 
  rounded <- round(pR2Sums, 3)
  text(bp, pR2Sums, labels = rounded, pos = 3, cex = 0.8)

}


par(mfrow = c(1,2))
pcpr2()
pcpr2(subset.instants = T)
type3mat
