# Brew targeted analysis

# Notes ----
# Data extracted from raw data files with ProFinder. A custom PCDL was used as a target for the extraction
#of 120 coffee compounds found in online databases. Once the data was extracted, each compound group was
#checked for poorly drawn peaks. Bad peaks were redrawn or removed if judged to be noise. Blanks were
#checked to make sure extracted peaks were not spurious. Any judged to be noise only were deleted.

# Selected samples, rejected samples and blanks were extracted (182 samples)
# Pooled QCs were extracted separately for calculation of CVs.

# 88 Compound groups were extracted from the coffee compounds initially collected. 17 were rejected
#for the following reasons, leaving 71.

#  1 due to detector saturation: caffeine
#  2 due to not of interest: Ochratoxin, 2 isomers
#  2 due to high single ion compounds: 2-Ac-5-Methylthiophene, 2-Me-5-(1-propenyl)pyrazine
#  8 due to high levels in blanks: Kahweol, 3-(3,4-dihydroxyphenyl)-2-propenoic acid, Linoleic acid, 
#        4-Ethyl-1,2-dimethoxybenzene, Campestanol, Trimethylamine, Feruloylquinic acid (3)
#  1 due to too many missing peaks: Liberine
#  1 due to mistaken identity: Trigonelline (RT is around 0.7 min)
#  3 due to high CV in pooled QCs: diacetin, Nicotinic acid, Phenylalanine

# Update at June 2018: More compounds excluded

# File was exported from ProFinder as "export detailed csv". Due to the commas in compound names, it was 
#difficult to parse the csv to be read in directly. Therefore, the csv was manually manipulated in Excel 
#to correct the shifted rows. Codes were given to each extracted compound eg cmpd_152_3.11.
#Excel file location: \\Inti\BMA\Coffee project NCI\Targeted analysis 88 cmpds Feb 15.xlsx.
#The "cleaned" sheet was exported to the R project directory.

# Data preparation 

# Readr automatically inserts NAs into missing cells.
# Read in and subset names and peak areas, subset selected coffees

# coffee metadata in Masked Batch List of coffee samples_data.files.xlsx


# Function to get scores and loadings ----
# Fig 1a publication: Targeted analysis PC1 vs PC2. First Log transform
# Fig 1b publication: Targeted analysis PC1 vs PC3

brew.pca <- function(filter = F, normalise = F) {
  
  library(tidyverse)
  library(pca3d)
  library(zoo)
  meta <- read_csv("Brew meta pt order march.csv")

  # Join intensities to metadata and filter selected samples
  raw  <- read_csv("Targeted 88 cmpds codes.csv") %>% left_join(meta, by = "data.file") %>% 
          select(plate:selected, data.file, everything()) %>% filter(selected == 1)
  
  # subset filter coffees or not
  if(filter == T) raw <- raw %>% filter(brew.method == "Filter") else raw

  # Read in compound codes
  codes <- read.csv("Compound codes.csv")

  # Filter rejected compounds (from above) and subset intensities only for PCA
  intsonly <- raw %>% select(-(plate:data.file), -ends_with("_rej")) %>% as.matrix
  
  # impute missings with median values
  mat <- na.aggregate(intsonly, FUN = median)
  
  # normalize to total intensity (normalised matrix needs transposing)
  normat <- apply(mat, 1, function(x) x/sum(x))
  logmat <- if(normalise == T) t(log2(normat)) else log2(mat)

  # Calculate PCs and specify coffee groups
  pcs <- prcomp(logmat, scale. = T)
  grps <- if(filter == T) factor(raw$roast) else factor(raw$brew.method)

  # get summmary of PCA for cumulative variance explained
  print(summary(pcs)) # First 3 PCs explain 80% (0.8033) variance

  # plot scores and loadings PC1 vs PC2, base R (not used, see old file version)
  # with pca3d for publication (first making data frame from pc object). Easier to show groups.
  # plots copied to clipboard at 941x638 px for draft
  plt <- pca2d(pcs, components = 1:2, group = grps, axe.titles = 
                         c("Score on PC1", "Score on PC2"))
  box(which = "plot", lty = "solid")
  legend("bottomright", plt$groups, col=plt$colors, pch=plt$pch)

  # Calculate contributions. Get absolute values of rotation and convert to proportions, join names, plot
  # do not use absolute values to get positive and negative contributions
  aload <- pcs$rotation
  contr <- sweep(aload, 2, colSums(aload), "/") %>% data.frame %>% 
                 rownames_to_column(var="cpd.code") %>% left_join(codes) %>%
                 select(cpd.code, cpd.name.new, everything())

  # plot top 10 contributions for PC1, 2, 3. First prepare data then plot
  contr.df <- if(filter == T) contr %>% select(cpd.name.new:PC1) else contr %>% select(cpd.name.new:PC3)
  df.plot <-  contr.df %>% 
              gather(PC, val, -cpd.name.new) %>% 
              mutate(absval = abs(val), sign = ifelse(val == absval, "Pos", "Neg")) %>% 
              group_by(PC) %>% top_n(n = 15, wt = absval)

  # Fill colour needs to go outside aes() to work
  base <-
    ggplot(df.plot, aes(x = val, y = cpd.name.new)) + 
    geom_segment(aes(x = 0, y = cpd.name.new, xend = val, yend = cpd.name.new), #size=1.5, 
                 lineend = "butt") + 
    theme_minimal() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point() + #ylab("Compounds contributing most to PC") + 
    xlab("Relative contribution to PC") +
    theme(legend.position = "none", plot.title = element_text(size = 18),
          axis.title.y = element_blank(),
          strip.text.x = element_text(hjust = 0, size = 10))
 
    brew.load <- if(filter == T) { 
      base + ggtitle("B")
        } else {
      base + facet_wrap( ~ PC, ncol = 1, scales = "free_y", shrink = T) + ggtitle("C")
        }
    
  return(brew.load)
}

# All data, unnormalised: scores and loadings plots
output <- brew.pca()
output

# Filter coffees only, normalised to remove effect of concentration
output.filt <- brew.pca(filter = T, normalise = T)
output.filt

# Correlation heatmap ----

brew.cor <- function() {
  meta <- read_csv("Brew meta pt order march.csv")

  #Join intensities to metadata and filter selected samples
  raw  <- read_csv("Targeted 88 cmpds codes.csv") %>% 
    left_join(meta, by = "data.file") %>% 
    select(plate:selected, data.file, everything()) %>% filter(selected == 1)

  #if(filter == T) raw <- raw %>% filter(brew.method == "Filter") else raw

  #Read in compound codes
  codes <- read.csv("Compound codes.csv")

  #Filter rejected compounds (from above) and subset intensities only for PCA
  intsonly <- raw %>% select(-(plate:data.file), -ends_with("_rej")) %>% as.matrix
  #Of all compounds (for annotation evidence)
  cormat <- cor(log2(intsonly), use = "pairwise.complete.obs")

  #make data frame of labels
  df <- data.frame(cpd.code = rownames(cormat)) %>% left_join(codes, by="cpd.code")
  
  #add rownames, remove colnames
  rownames(cormat) <- df$abrev.new
  colnames(cormat) <- rep(NA, ncol(intsonly))

  #rownames(cormat) <- paste(rownames(cormat), 1:71, sep = ", ")
  library(corrplot)
  cplot <- corrplot(cormat, method="square", tl.col="black", order="hclust", tl.cex=0.7,  
         cl.ratio=0.1, cl.align="r", cl.pos="r", cl.cex=0.8)
  return(cplot)

}
brew.cor()

# Total intensities by brew method ----

library(dplyr)
meta <- read_csv("Brew meta pt order march.csv")
table(meta$brew.method)

# Join intensities to metadata and filter selected samples
raw <- read_csv("Targeted 88 cmpds codes.csv") %>% left_join(meta, by = "data.file") %>% 
  filter(selected == 1) %>% select(brew.method:roast, contains("cpd_"), -contains("_rej"))

df <- raw %>% transmute(total.int = (rowSums(select(., -1), na.rm = T)/1000000)) %>% 
  bind_cols(raw[, 1])

library(gplots)
boxplot2(df$total.int ~ df$brew.method, varwidth = T, col= "limegreen", 
         top =T, ylab="Total intensity (millions of counts)", main = 
           "Targeted analysis (64 compounds)")

# t-test to test difference between espressos and instants
t.test(df$total.int ~ df$brew.method, subset = df$brew.method == "Instant" | 
         df$brew.method == "Espresso")
