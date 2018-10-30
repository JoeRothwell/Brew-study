#Read in heatmap data. .csv made from Data for coffee heatmap.xlsx sheet
#Coffees were ordered by characteristic for easier comparison

coffees <- read.csv("Brew targeted 93 cmpds.csv", header=T)
mat     <- as.matrix(coffees[, 2:94])
rownames(mat) <- coffees[, 1]

#replace missing values (1s) with half the minimum value for each cmpd (row)
halfmin <- function(x) replace(x, x == 1, min(x[x != 1])/2)
mat2    <- apply(mat, 2, halfmin)

#log transform for normal distribution
logmat <- log2(mat2)

#scan in character vectors of labels from csv files. Names that start with numbers need to
#be surrounded by "" in the csv for scan to work
cmpdnames <- scan("coffee compound labels.csv", what="character", sep=",")

library(gplots)
heatmap.2(t(logmat), 
  Rowv = T, Colv = T, col = redblue(256), 
  trace = "none", 
  dendrogram = "both",
  density.info="none",
  keysize = 1,
  scale = "row", key = T,
  labRow = cmpdnames,
  margins = c(8,11), 
  offsetCol = 0.005, offsetRow = 0.005, 
  cexRow = 0.4, cexCol = 0.4
  #coffee numbers are 18 other 27 filter 6 other 12 instant 13 other
  #ColSideColors = c(rep("gray", 18), rep("green", 27), rep("gray", 6), 
                            #rep("yellow", 12), rep("gray", 13)),
  )

length(logmat)
hist(logmat, breaks=100, col=3)

library(corrplot)
cmat <- cor(logmat, method = "pearson")
rownames(cmat) <- cmpdnames
corrplot(cmat, method="square", tl.col="black", order="hclust", 
         tl.cex=0.5, cl.ratio=0.1, cl.align="r",
         cl.pos="r", cl.cex=0.5, mar=c(1,1,2,1))

#pca
pca <- prcomp(logmat, scale. = T)
plot(pca$x[, 1], pca$x[, 2], pch=19, col=2, cex=0.5)
biplot(pca, cex=0.5)
