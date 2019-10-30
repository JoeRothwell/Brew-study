# Exploratory analysis of coffee brews. Metadata needs to be in feature table order
# \\inti\BMA\Coffee Project\NCI\Sample lists\Masked Batch List of coffee samples_data files.xlsx

library(tidyverse)
cofints <- read_csv("Brew PT Feb 2015.csv", skip=8)
brews   <- read_csv("Brew meta pt order march.csv")

# Function to prep data and do PCA with log and/or scaling. Also returns table of PC scores.
brew.pca <- function(ft, meta, logtr = T, ...) {
 
  library(dplyr) 
  library(MetabolAnalyze)

  ints <- ft %>% select(ends_with("(raw)")) %>% t
  meta$replicate <- as.factor(meta$replicate)
  
  #write.csv(ints, "Peak table brew study.csv")
  
  # Subset metadata and intensities to be used (remove blanks for PCA)
  samples <- meta$sample.type != "Blank"
  mat <- ints[samples, ]
  
  #replace matrix zeros or NAs with 1s
  mat <- ifelse(mat == 0 | is.na(mat), 1, mat)
  
  #remove zero variance columns
  ifelse(apply(mat, 2, var) == 0, mat <- mat[, apply(mat, 2, var) != 0], mat)
  
  #log transform, scale, PCA
  mat <- if(logtr == T) log2(mat) else mat
  scalemat <- scaling(mat, ...)
  pc <- prcomp(scalemat, scale. = F, rank. = 10)
  
  df <- data.frame(meta[ samples, ], pc$x) #%>%
    #select(replicate, PC1:PC3) %>% arrange(replicate)
  
}
output <- brew.pca(cofints, brews, logtr = T, type = "pareto")

ggplot(output, aes(x=PC1, y=PC2, colour=brew.method)) + #geom_point() + 
  geom_text(aes(label = replicate)) + 
  theme_bw() + 
  geom_hline(yintercept = 0, colour="grey") + 
  geom_vline(xintercept = 0, colour="grey")


# Dendrogram (to replace with circular dendrogram)
mat <- cofints %>% select(ends_with("(raw)")) %>% t %>% log2
  
#log transform and scale
scalemat <- scaling(mat, type = "pareto")
hh1 <- hclust(dist(scalemat))
brew_method <- as.numeric(factor(brews$brew.method))
ColorDendrogram(hh1, y=brew_method, branchlength=38, labels=brews$brew.method)


# Total usable signal ----

# Get total usable signal of the selected coffees
ints <- ft %>% select(ends_with("(raw)")) %>% t
ints.filt <- ints[meta$selected == 1, ]
totalints <- rowSums(ints.filt)
totalints <- data.frame(meta[meta$selected == 1, ], 
             totalint = (as.numeric(rowSums(ints.filt)/1000000) ))

# Plot boxplot of intensities by coffee brew. Set margins and plot for 2 figures
par(mar = c(2,2,3,2)) 
par(mfrow = c(2,1))

library(gplots)
boxplot2(totalints$totalint ~ totalints$brew.method, main = "Untargeted analysis (all features)",
        ylab = "Total usable intensity (millions of counts)", varwidth = T, col = "dodgerblue", top = T)

t.test(totalints$totalint ~ totalints$brew.method, 
       subset = totalints$brew.method == "Instant" | totalints$brew.method == "Espresso")

# Feature counts ----
brew <- bind_cols(meta, data.frame(ints))

# find number of features present in X coffees
brew.rep1 <- brew %>% filter(replicate == 1) %>% select(X1:ncol(brew))
nfeatures <- colSums(brew.rep1 > 1)

hist(nfeatures, xlab="No. coffees", ylab="No. features", breaks=100, col="lightblue")
plot(ecdf(nfeatures), xlab="No. coffees", ylab="No. features")

# Volcano plots----
# The following code is old and may not work
# Normal v decaf and Arabica vs blend
# get logical vectors for the 4 groups under study. Check with sum(normal==T)
normal  <- meta$caffeine  == "Caf"     & meta$selected == T #n=60
decafs  <- meta$caffeine  == "Decaf"   & meta$selected == T #n=16
arabica <- meta$bean.type == "Arabica" & meta$selected == T #n=38
blend   <- meta$bean.type == "Blend"   & meta$selected == T #n=31

# normal/decaf fold changes and t-test
fc.caffeine <- apply(ints[normal, ], 2, mean) / apply(ints[decafs, ], 2, mean)
tt.caffeine <- apply(log2(ints), 2, function(x) t.test(x[normal], x[decafs]))
praw <- sapply(tt.caffeine, "[", c(3)) %>% unlist
padj <- p.adjust(praw, method="BH")
plot(log2(fc.caffeine), log10(padj), ylim=c(0, -40))

# Arabica/blend fold changes and t-test
fc.beantype <- apply(ints[arabica, ], 2, mean) / apply(ints[blend, ], 2, mean)
tt.beantype <- apply(log2(ints), 2, function(x) t.test(x[arabica], x[blend]))
praw2 <- sapply(tt.caffeine, "[", c(3)) %>% unlist
padj2 <- p.adjust(p, method="BH")
plot(log2(fc.beantype), log10(padj), ylim=c(0, -40))

# Could also plot feature boxplots eg
cafeffect <- (meta$caffeine == "Caf" | meta$caffeine == "Decaf") & meta$selected == T


