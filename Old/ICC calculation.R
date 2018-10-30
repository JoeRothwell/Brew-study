# 03/11/2016
# Mayo Fecal Metabolomics Methods Study - estimating precision, accuracy, and stability of fecal metabolites
# Programer: Erikka Loftfield

library(lme4)
?lmer

#### Get working directory
getwd()
setwd("H:/Mayo/SAS/output/data/LOD/new")

####What if I want the "precision"

ICC.fun1 <- function(data, id, X) {
  id2    <- data[, id]
  tm     <- lmer(data[, X] ~ 1 + (1|id2))
  s2e    <- (attr(VarCorr(tm), "sc")) ^2
  s2sub  <- VarCorr(tm)$id2[1]
  s2tot  <- s2e + s2sub
  ICC    <- s2sub/s2tot
  ICC
}

#####Fresh


#### Read a .csv file with all data
fresh_lod_log <- read.table(file="fresh_lod_log.CSV", sep=",", header=TRUE)
fresh_lod     <- read.table(file="fresh_lod.CSV", sep=",", header=TRUE)

head(fresh_lod, n=2)
names(fresh_lod)


#nc hold only met variables
nc_fresh_log <- ncol(fresh_lod_log) - 2
nc_fresh     <- ncol(fresh_lod) - 2

#fresh_log_ <-log(fresh_lod[c(1:625)])

####calculate ICCs with all of the data --- FOR MAIN RESULTS
fresh_lod_log$ID        <- as.character(fresh_lod_log$ID)
fresh_lod_log$CONDITION <- as.character(fresh_lod_log$CONDITION)

fresh_lod$ID        <- as.character(fresh_lod$ID)
fresh_lod$CONDITION <- as.character(fresh_lod$CONDITION)

ICCVec_fresh <- rep(NA, nc_fresh)

for(i in 1:nc_fresh){
  ICCVec_fresh[i] <- ICC.fun1(fresh_lod, "ID", colnames(fresh_lod)[i])
  cat(i, colnames(fresh_lod)[i], ICCVec_fresh[i], "\n")
}

hist(ICCVec_fresh, 100, main="Precision of Fresh Frozen (D0) Samples", xlab="ICC", plot=TRUE) 
summary(ICCVec_fresh)
