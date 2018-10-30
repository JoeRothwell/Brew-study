#PCA and exploratory analysis (this script to be put in correct order)
#Original coffee metadata was imported from:
#\\inti\BMA\Coffee Project\NCI\Sample lists\Masked Batch List of coffee samples_data files.xlsx
#put in peak table order for match with intensity data
library(tidyverse)
library(MetabolAnalyze)

# Exploratory analysis ---------------------------------------------------------------------------

#read in brew intensity data and metadata (200 obs incl blanks and QCs)
ints <- read_csv("Brew PT Feb 2015.csv", skip=8) %>% select(2:201) %>% t
#write.csv(ints, "Peak table brew study.csv")

#meta <- read_csv("Brew labels pt order feb plate.csv")
meta <- read_csv("Brew meta pt order march.csv")

#define function to prep data and do PCA with log and/or scaling. Also returns table of PC scores.
brew.explore <- function(samples=c(2:201), logtr = TRUE, ...) {
  
  #exploratory analysis. First define metadata and intensities to be used
  mat <- ints[samples, ]
  
  #replace matrix zeros or NAs with 1s
  mat <- ifelse(mat == 0 | is.na(mat), 1, mat)
  
  #remove zero variance columns
  ifelse(apply(mat, 2, var) == 0, mat <- mat[, apply(mat, 2, var) != 0], mat)
  
  #log transform and scale
  if(logtr == T) logmat <- log2(mat) else logmat <- mat 
  scalemat <- scaling(mat, ...)
  
  #calculate PCA and plot scores
  pc <- prcomp(scalemat, scale. = F)
  #plot(pc$x[, 1], pc$x[, 2], xlab = "Scores on PC1", ylab = "Scores on PC2")
  #abline(v=0, h=0, lty=3)
  
  df <- data.frame(meta[ samples, ], pc$x)
  p <- ggplot(df[, samples], aes(x=PC1, y=PC2, colour=replicate.code, shape=caffeine)) + 
    geom_text(aes(label = replicate)) + 
    #geom_point() + 
    theme_bw() + geom_hline(yintercept = 0, colour="grey") + 
    geom_vline(xintercept = 0, colour="grey")
  print(p)
  
  df2 <- df %>% select(replicate, PC1:PC3) %>% arrange(replicate)
  return(data.frame(df2))
  
}

#Run function removing blanks at run index 67, etc. Scaling types are "unit", "pareto", "none". 
result <- brew.explore(samples=c(-67,-76,-(177:182)), logtr = T, type = "unit")

brew.explore(logtr = T, type = "unit")

#Define function to perform hierarchical clustering
brew.cluster <- function(logtr = TRUE, ...) {
  
  #exploratory analysis. First define metadata and intensities to be used
  mat <- t(ints) 
  
  #log transform and scale
  if(logtr == T) mat <- log2(mat) else mat <- mat 
  scalemat <- scaling(mat, ...)
  
  hh1 <- hclust(dist(scalemat))
  caffeine_as_numeric = as.numeric(factor(meta$brew.method))
  ColorDendrogram(hh1, y=caffeine_as_numeric, branchlength=38, labels=meta$brew.method)
  
}
brew.cluster(type="unit")

# Total usable signal --------------------------------------------------------------------------

#Get total usable signal of the selected coffees
ints.filt <- ints[meta$selected == 1, ]
totalints <- rowSums(ints.filt)
totalints <- data.frame(meta[selected == 1, ], totalint = (as.numeric(rowSums(ints.filt)/1000000) ))

#Plot boxplot of intensities by coffee brew. Set margins and plot for 2 figures
par(mar=c(2, 2, 3, 2)) 
par(mfrow=c(2,1))

library(gplots)
boxplot2(totalints$totalint ~ totalints$brew.method, main = "Untargeted analysis (all features)",
        ylab = "Total usable intensity (millions of counts)", varwidth = T, col = "dodgerblue", top = T)

t.test(totalints$totalint ~ totalints$brew.method, 
       subset = totalints$brew.method == "Instant" | totalints$brew.method == "Espresso")



# Feature counts ---------------------------------------------------------------------------------
brew <- bind_cols(meta, data.frame(ints))

#find number of features present in X coffees
brew.rep1 <- brew %>% filter(replicate == 1) %>% select(X1:ncol(brew))
nfeatures <- colSums(brew.rep1 > 1)

hist(nfeatures, xlab="No. coffees", ylab="No. features", breaks=100, col="lightblue")
plot(ecdf(nfeatures), xlab="No. coffees", ylab="No. features")

#Top features: read in intensity data and metadata

#subset intensity data for selected observations
mat  <- ints[ which(meta$selected == 1), ]
meta <- meta[meta$selected == 1, ]

#loop through each feature getting the highest median
topfeat <- function(cofvec){
  #make an empty vector of length 10
  output <- numeric(10)
  #copy matrix for resetting later
  resetmat <- mat
  #find the index of the feature with highest median value and assign to top
  for (i in 1:10) {
    top              <- which.max(apply(mat[cofvec, ], 2, median))
    output[i]        <- top     #assign this value to the empty vector
    mat[cofvec, top] <- 0       #assign zero to this feature so the next highest can be found with max
  }
  mat <- resetmat  #reset the matrix
  return(output)
}

#run function on coffee subsets
all      <- topfeat(c(1:76))
normal   <- topfeat(meta$caffeine=="Caf")
decaf    <- topfeat(meta$caffeine=="Decaf")
arabica  <- topfeat(meta$bean.type=="Arabica")
blend    <- topfeat(meta$bean.type=="Blend")
boiled   <- topfeat(meta$brew.method=="Boiled")
coldbrew <- topfeat(meta$brew.method=="Cold Brew")
espresso <- topfeat(meta$brew.method=="Espresso")
frpress  <- topfeat(meta$brew.method=="French Press")
instant  <- topfeat(meta$brew.method=="Instant")
kcup     <- topfeat(meta$brew.method=="K-Cup")
percol   <- topfeat(meta$brew.method=="Percolated")
light    <- topfeat(meta$roast=="Light")
med      <- topfeat(meta$roast=="Med")
dark     <- topfeat(meta$roast=="Dark")
allcofs  <- cbind(all, normal, decaf, arabica, blend, boiled, coldbrew, espresso, frpress, 
                  instant, kcup, percol, light, med, dark)

#Alternative: using a loop in a loop to subset coffee groups. Loop through coffee groups and features
#First get ranks for coffee attributes of interest
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
      top <-
        which.max(apply(mat[meta[, attribute] == j,] , 2, median))
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

#---------------------------------------------------------------------------------------------

#volcano plots, Normal v decaf and Arabica vs blend
#get logical vectors for the 4 groups under study. Check with sum(normal==T)
normal  <- meta$caffeine  == "Caf"     & meta$selected == T #n=60
decafs  <- meta$caffeine  == "Decaf"   & meta$selected == T #n=16
arabica <- meta$bean.type == "Arabica" & meta$selected == T #n=38
blend   <- meta$bean.type == "Blend"   & meta$selected == T #n=31

#normal/decaf fold changes and t-test
fc.caffeine <- apply(ints[normal, ], 2, mean) / apply(ints[decafs, ], 2, mean)
tt.caffeine <- apply(log2(ints), 2, function(x) t.test(x[normal], x[decafs]))
praw <- sapply(tt.caffeine, "[", c(3)) %>% unlist
padj <- p.adjust(praw, method="BH")
plot(log2(fc.caffeine), log10(padj), ylim=c(0, -40))

#Arabica/blend fold changes and t-test
fc.beantype <- apply(ints[arabica, ], 2, mean) / apply(ints[blend, ], 2, mean)
tt.beantype <- apply(log2(ints), 2, function(x) t.test(x[arabica], x[blend]))
praw2 <- sapply(tt.caffeine, "[", c(3)) %>% unlist
padj2 <- p.adjust(p, method="BH")
plot(log2(fc.beantype), log10(padj), ylim=c(0, -40))

#could also plot feature boxplots eg
cafeffect <- (meta$caffeine == "Caf" | meta$caffeine == "Decaf") & meta$selected == T


