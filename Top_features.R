# Read in intensity data and metadata
ints <- read.csv("Brew PT Feb 2015.csv", skip = 8, row.names=1)
meta <- read.csv("Brew meta pt order march.csv")

# Subset intensity data for selected observations
mat  <- t(ints[, which(meta$selected == 1)])
meta <- meta[meta$selected == 1, ]

# Loop through each feature getting the highest median
topfeat <- function(cofvec){
  
  #make an empty vector of length 10
  output   <- numeric(10)
  resetmat <- mat #copy matrix for resetting later
  
  for (i in 1:10) { 
    #find the index of the feature with highest median value and assign to top
    top <- which.max(apply(mat[cofvec, ], 2, function(x) median(x, na.rm = T)))
    output[i] <- top       #assign this value to the empty vector
    mat[cofvec, top] <- 0  #assign zero to this feature so the next highest can be found with max
  }
  
mat <- resetmat  #reset the matrix
return(output)
}

all      <- topfeat(c(1:76))
normal   <- topfeat(meta$caffeine == "Caf")
decaf    <- topfeat(meta$caffeine == "Decaf")
arabica  <- topfeat(meta$bean.type == "Arabica")
blend    <- topfeat(meta$bean.type == "Blend")
boiled   <- topfeat(meta$brew.method == "Boiled")
coldbrew <- topfeat(meta$brew.method == "Cold Brew")
espresso <- topfeat(meta$brew.method == "Espresso")
frpress  <- topfeat(meta$brew.method == "French Press")
instant  <- topfeat(meta$brew.method == "Instant")
kcup     <- topfeat(meta$brew.method == "K-Cup")
percol   <- topfeat(meta$brew.method == "Percolated")
light    <- topfeat(meta$roast == "Light")
med      <- topfeat(meta$roast == "Medium")
dark     <- topfeat(meta$roast == "Dark")
allcofs  <- cbind(all, normal, decaf, arabica, blend, boiled, coldbrew, espresso, frpress, 
             instant, kcup, percol, light, med, dark)

# Alternative:---- 
# A loop in a loop to subset coffee groups. Loop through coffee groups and features
# First get ranks for coffee attributes of interest
library(dplyr)
meta <- meta %>% filter(selected == 1) %>% mutate(CafID = dense_rank(caffeine), 
            btID = dense_rank(bean.type), bmID = dense_rank(brew.method), allID = 1)

top.feat <- function(meta, attribute) {
  outputall <-
    vector(mode = "list", length = length(unique(meta[, attribute])))
  for (j in unique(meta[, attribute])) {
    #outer loop (cycles through subsets)
    output <- numeric(10)
    resetmat <- mat
    for (i in 1:10) {
      #inner loop (gets top features)
      top <- which.max(apply(mat[meta[, attribute] == j,] , 2, 
              function(x) median(x, na.rm = T)))
      output[i] <- top
      names(output)[i] <- colnames(mat)[top]
      mat[(meta[, attribute] == j), top] <- 0
    }
    mat <- resetmat #matrix is reset after each subset
    outputall[[j]] <- output
    names(outputall)[j] <- paste(attribute, j, sep=".")
  }
  return(outputall)
}

alllist <- top.feat(meta, attribute="allID")
caflist <- top.feat(meta, attribute="CafID")
btlist  <- top.feat(meta, attribute="btID")
bmlist  <- top.feat(meta, attribute="bmID")

rbind(unlist(caflist), unlist(btlist), unlist(bmlist))
