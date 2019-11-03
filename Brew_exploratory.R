# Exploratory analysis of coffee brews. Metadata needs to be in feature table order
# \\inti\BMA\Coffee Project\NCI\Sample lists\Masked Batch List of coffee samples_data files.xlsx

library(tidyverse)
ints <- read_csv("Brew PT Feb 2015.csv", skip=8)
brews <- read_csv("Brew meta pt order march.csv")

# Function to prep data and do PCA with log and/or scaling. Also returns table of PC scores.
brew.pca <- function(ft, meta, logtr = T, ...) {

  ints <- ft %>% select(ends_with("(raw)")) %>% t
  meta$replicate <- as.factor(meta$replicate)
  
  # Subset metadata and intensities to be used (remove blanks for PCA)
  samples <- meta$sample.type != "Blank"
  mat <- ints[samples, ]
  
  # replace matrix zeros or NAs with 1s, remove zero variance cols
  mat <- ifelse(mat == 0 | is.na(mat), 1, mat)
  ifelse(apply(mat, 2, var) == 0, mat <- mat[, apply(mat, 2, var) != 0], mat)
  
  #log transform, scale, PCA
  mat <- if(logtr == T) log2(mat) else mat
  scalemat <- MetabolAnalyze::scaling(mat, ...)
  pc <- prcomp(scalemat, scale. = F, rank. = 10)
  
  df <- data.frame(meta[ samples, ], pc$x)

}
output <- brew.pca(ints, brews, logtr = T, type = "pareto")

ggplot(output, aes(x=PC1, y=PC2, colour=brew.method)) + #geom_point() + 
  geom_text(aes(label = replicate)) + 
  theme_bw() + 
  geom_hline(yintercept = 0, colour="grey") + 
  geom_vline(xintercept = 0, colour="grey")

# Dendrogram (to replace with circular dendrogram)
mat <- ints %>% select(ends_with("(raw)")) %>% t %>% log2
  
#log transform and scale
library(MetabolAnalyze)
scalemat <- scaling(mat, type = "pareto")
hh1 <- hclust(dist(scalemat))
brew_method <- as.numeric(factor(brews$brew.method))
library(sparcl)
ColorDendrogram(hh1, y=brew_method, branchlength=38, labels=brews$brew.method)

library(dendextend)
dend1 <- as.dendrogram(hh1) %>%
  set("branches_k_color", k=4) %>% 
  #set("labels_col", cmpds) %>%
  set("labels_cex", c(0.6)) #%>%
#hang.dendrogram(hang_height = 0.3)

circlize_dendrogram(dend1, dend_track_height = 0.85)

# Total usable signal ----

# Get total usable signal of the selected coffees
ints <- ints %>% select(ends_with("(raw)")) %>% t
ints1 <- ints[brews$selected == 1, ]
totals <- rowSums(ints1)
totals <- data.frame(brews[brews$selected == 1, ], 
             totalint = (as.numeric(rowSums(ints1)/1000000) ))

# Plot boxplot of intensities by coffee brew. Set margins and plot for 2 figures
par(mar = c(2,2,3,2)) 
par(mfrow = c(2,1))

library(gplots)
boxplot2(totals$totalint ~ totals$brew.method, main = "Untargeted analysis (all features)",
        ylab = "Total usable intensity (millions of counts)", varwidth = T, 
        col = "dodgerblue", top = T)

t.test(totals$totalint ~ totals$brew.method, 
       subset = totals$brew.method == "Instant" | totals$brew.method == "Espresso")

# Feature counts ----
brew <- bind_cols(brews, data.frame(ints))

# find number of features present in X coffees
rep1 <- brew %>% filter(replicate == 1) %>% select(X1:ncol(brew))
nfeatures <- colSums(rep1 > 1)

hist(nfeatures, xlab="No. coffees", ylab="No. features", breaks=100, col="lightblue")
plot(ecdf(nfeatures), xlab="No. coffees", ylab="No. features")

# Volcano plots (old) ----
# Normal v decaf and Arabica vs blend
# get logical vectors for the 4 groups under study. Check with sum(normal)
meta <- brews
normal  <- meta$caffeine  == "Caf"     & meta$selected == T  #n=60
decafs  <- meta$caffeine  == "Decaf"   & meta$selected == T   #n=16
arabica <- meta$bean.type == "Arabica" & meta$selected == T & !is.na(meta$bean.type) #n=38
blend   <- meta$bean.type == "Blend"   & meta$selected == T & !is.na(meta$bean.type) #n=31

# normal/decaf fold changes and t-test
fc1 <- apply(ints[normal, ], 2, mean) / apply(ints[decafs, ], 2, mean)
tt1 <- apply(log2(ints), 2, function(x) t.test(x[normal], x[decafs]))
praw <- sapply(tt1, "[", c(3)) %>% unlist
padj <- p.adjust(praw, method="BH")
plot(log2(fc1), log10(padj), ylim=c(0, -40))

# Arabica/blend fold changes and t-test
fc2 <- apply(ints[arabica, ], 2, mean) / apply(ints[blend, ], 2, mean)
tt2 <- apply(log2(ints), 2, function(x) t.test(x[arabica], x[blend]))
praw2 <- sapply(tt2, "[", c(3)) %>% unlist
padj2 <- p.adjust(praw2, method="BH")
plot(log2(fc2), log10(padj), ylim=c(0, -40))

# Could also plot feature boxplots eg
cafeffect <- (meta$caffeine == "Caf" | meta$caffeine == "Decaf") & meta$selected == T


# PCPR2 ----
# on either all coffee samples (n=174) or all coffee samples except instants (n=150)
# Final analysis for coffee brew paper
# Run PCPR2_coffee.R lines 7-38

# Load the data in the X matrix containing NMR spectra and Z matrix containing the list of explanatory variables of interest
ints <- read_csv("Brew PT Feb 2015.csv", skip = 8) %>% select(2:201) %>% t
#"Unknown" was replaced with NA in metadata to avoid aliased coefficients in model

#174 coffees to be included (exclusion of blanks and QCs) or 150 (exclusion of instants)
#samples <- !is.na(meta$coffee.brew)
samples <- !is.na(meta$coffee.brew) & meta$brew.method != "Instant"

#features detected in more than 2 samples only
X_DataMatrix <- ints[samples, ]
filt <- colSums(X_DataMatrix > 1) > 2

#subset both dimensions. Leaves 3626 and 3572 features for full and subset data.
X_DataMatrix <- log2(ints[samples, filt])
Z_InterestFactors <- meta[samples, c(1:2, 10:13)]
X_DataMatrixCentered = scale(X_DataMatrix, center = TRUE, scale = FALSE)
X_DataMatrixScaled = scaling(X_DataMatrixCentered, type = "pareto")

library(pcpr2)
output <- runPCPR2(X_DataMatrixScaled, Z_InterestFactors)
plot(output)
summary(output)


