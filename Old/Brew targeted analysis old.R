#Brew targeted analysis

#Data extracted from raw data files with ProFinder. A custom PCDL was used as a target for the extraction
#of 120 coffee compounds found in online databases. Once the data was extracted, each compound group was
#checked for poorly drawn peaks. Bad peaks were redrawn or removed if judged to be noise. Blanks were
#checked to make sure extracted peaks were not spurious. Any judged to be noise only were deleted.

#Selected samples, rejected samples and blanks were extracted (182 samples)
#Pooled QCs were extracted separately for calculation of CVs.

#88 Compound groups were extracted from the coffee compounds initially collected.
#Rejected due to high single ion compounds: 2-Ac-5-Meth, 2-Me-5-(1-propenyl)pyrazine 
#Rejected due to high levels in blanks: Kahweol, 3-(3,4-dihydroxyphenyl)-2-propenoic acid, Linoleic acid, 
#4-Ethyl-1,2-dimethoxybenzene, Campestanol, Trimethylamine, Feruloylquinic acid (3) (deleted from XL sheet)
#Rejected due to saturation: caffeine

#File was exported from ProFinder as "export detailed csv". Due to the commas in compound names, it was 
#difficult to parse the csv to be read in directly. Therefore, the csv was manually manipulated in Excel 
#to correct the shifted rows. Also, some compound names were simplified at the same time. The Excel file 
#is stored #in \\Inti\BMA\Coffee project NCI\Targeted analysis 88 cmpds Feb 15.xlsx.
#The "cleaned" sheet was exported to the R project directory.

#Readr automatically inserts NAs into missing cells.
#Read in and subset names and peak areas, subset selected coffees

#coffee metadata in Masked Batch List of coffee samples_data.files.xlsx
library(tidyverse)
meta <- read_csv("Brew meta pt order feb.csv")
raw  <- read_csv("Targeted 81 cmpds Mar 15.csv") %>% left_join(meta, by = "data.file") %>% 
  select(plate:selected, data.file, everything())

#subset selected samples and blanks or filter coffees only
sel <- raw %>% filter(selected == 1) # or
sel <- raw %>% filter(selected == 1 & brew.method == "Filter")
samp <- raw %>% filter(sample.type == "Sample" | sample.type == "Defat_sample")
bls <- raw %>% filter(sample.type == "Blank")
filt <- raw %>% filter(selected == 1 & brew.method == "Filter")

#----------------------------------------------------------------------------------------------
#Compound subsets
#1. Determine highest intensity compounds from selected replicates. First subset intensities only
intonly  <- sel %>% select(-(plate:data.file))
cmpdvec2 <- rank(colSums(intonly, na.rm = T)) > 60
TOPcmpds <- intonly[ , cmpdvec2]

#2. Determine NCI compounds of interest (see document from Erikka)
#Trigonelline, nicotinic acid, caffeine + 3 metabolites, caffeic acid, caffeoylquinic acid (3 isomers), 
#dicaffeoylquinic acid (2 isomers), ferulic acid, feruloylquinic acid (3 isomers), coumaric acid, 
#coumaroylquinic acid (2 isomers), cafestol, kahweol, linoleic acid, hydroxybenzoic acid (2 isomers), 
#hydroxyglutarate, phenylalanine, cinnamoylglycine, cyclo(leu-pro).

#cmpdvec1 <- c(33,46,47,49:53,60,62:65,74,75,77:80,83,84,86,88,92,95,96,98:100)
cmpdvec1 <- c(32,44:50,56,58:61,70,71,73:76,79,80,82,86,89,90,92:94)
NCIcmpds <- sel[ , cmpdvec1]

#3. Subset filter coffees only
intonly.filt <- filt %>% select(-(plate:data.file))

#Melt data and remove NAs, plot boxplot
ints  <- NCIcmpds %>% gather(Compound, intensity) %>% na.omit
ints2 <- TOPcmpds %>% gather(Compound, intensity) %>% na.omit

ggplot(ints2, aes(y = log10(intensity), x = reorder(factor(Compound), intensity, median))) + 
  geom_boxplot(fill = "grey") + coord_flip() + theme_bw() +
  xlab("Compound") + ggtitle("NCI compounds of interest") #+ ggtitle("Most intense compounds")

#--------------------------------------------------------------------------------------------

#PCA
#Log transform and impute missings with median values
logmat <- log(intonly.filt)
logmatimp <- apply(logmat, 2, function(x){
  x[which(is.na(x))] <- median(x, na.rm=T)
  return(x)
})

#Calculate PCs and plot
pcs <- prcomp(logmatimp, scale. = T)
grps <- factor(sel$brew.method)
plot(pcs$x[, 1], pcs$x[, 2], col = grps, pch=19)
plot(pcs$rotation[, 1], pcs$rotation[, 2])
#legend
legend(7, 4.3, unique(grps), col = 1:length(grps), pch=19)

#For filter coffees only
pcs <- prcomp(logmatimp, scale. = T)
grps <- factor(filt$roast)

#Calculate contributions. Get absolute values of rotation and convert to proportions, then plot
pcs$rotation
aload <- abs(pcs$rotation)
contributions <- sweep(aload, 2, colSums(aload), "/")
par(mai=c(1,2,1,1))
barplot(contributions[, 1], horiz = T, las = 1)
ggplot(contributions, aes(y = PC1)) + geom_bar()

library(pca3d)
pca2d(pcs, components = 1:2, group = grps, legend = "bottom")
plot(pcs)
plot(pcs$rotation[, 1])

#heatmap
library(gplots)
heatmap.2(t(logmatimp), 
          Rowv = T, Colv = T, col = redblue(256), 
          trace = "none", 
          dendrogram = "both",
          density.info="none",
          keysize = 1,
          scale = "row", key = T,
          #labRow = cmpdnames,
          margins = c(2,9), 
          offsetCol = 0.005, offsetRow = 0.005, 
          cexRow = 0.7, cexCol = 0.4
          #27 filters, 12 instants, 9 espressos, 6 cold brews
          #ColSideColors = c(rep("gray", 18), rep("green", 27), 
          #rep("gray", 6), rep("yellow", 12), rep("gray", 13))
)

#correlation heatmap
cormat <- cor(logmat, use = "pairwise.complete.obs")
colnames(cormat) <- NULL
library(corrplot)
corrplot(cormat, method="square", tl.col="black", order="hclust", tl.cex=0.7, 
         cl.ratio=0.1, cl.align="r", cl.pos="r", cl.cex=0.8)

library(gplots)
heatmap.2(logmatimp, trace = "none")

#------------------------------------------------------------------------------------------------
#boxplot of individual features
boxplot(sel$Cafestol ~ sel$brew.method, varwidth=T, col="grey")
boxplot(sel$Scopoletin ~ sel$brew.method, col="grey")
boxplot(sel$'Benzenetriol (ii)' ~ sel$brew.method, col="grey")
boxplot(sel$Trigonelline ~ sel$brew.method, col="grey")

#icc: use samp subset
library(lme4)
lmer(Cafestol ~ brew.method + 1|coffee.brew, data = samp)
