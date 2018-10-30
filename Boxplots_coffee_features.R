#plot data on coffee types. read in and summarise data
brewlabs <- read.csv("Brew labels peak table order.csv")
#subset one replicate of each coffee and omit QCs and blanks
cf <- brewlabs[brewlabs$duplicate.no == 1, ]

#plot to get overview of data
library(ggplot2)
library(reshape2)
cftypes <- cf[, c(14:17)]
#dcast(cftypes, caffeine)
ggplot(cf, aes(x=caffeine, fill=brew.method)) + geom_bar(colour="black", position="dodge") + theme_grey() 
ggplot(cf, aes(x=roast, fill=caffeine)) + geom_bar(colour="black") + theme_bw()
ggplot(cf, aes(x=bean.type, fill=caffeine)) + geom_bar(colour="black") + theme_bw()
ggplot(cf, aes(x=brew.method, fill=caffeine)) + geom_bar(colour="black") + theme_bw()

brew <- read.csv("Joe Brew.csv", skip=8)
mat <- data.matrix(brew[ , (2:201)])
rownames(mat) <- brew$Mass
cfdata <- mat[ , brewlabs$duplicate.no == 1]
cfn <- droplevels(cf$caffeine)
roast <- droplevels(cf$roast)

#histogram distribution of feature count over 74 coffees
hist(apply(cfdata, 2, function(x) sum(x > 1)))


#distribution of numbers of coffees in which x features are present
hist(apply(cfdata, 1, function(x) sum(x > 1)), breaks=30) #most features are present in 60 coffees or more
plot(ecdf(apply(cfdata, 1, function(x) sum(x > 1)))) # ecdf
sum(apply(cfdata, 1, function(x) sum(x > 1) == 74)) # 1137 features in every coffee

#features unique to certain coffee classes: first subset data from logical vectors
cafs <- cfdata[ , cf$caffeine == "Caf" ]
decafs <- cfdata[ , cf$caffeine == "Decaf" ]

#find features which are absent from both subsets separately
cafzeros <- which(apply(cafs, 1, function(x) sum(x > 1) == 0))
decafzeros <- which(apply(decafs, 1, function(x) sum(x > 1) == 0))
zerocafions <- setdiff(cafzeros, decafzeros)
zerodecafions <- setdiff(decafzeros, cafzeros)

#plot boxplots of features. Convert first to dataframe and bind to caffeination factor, then melt
cdf1 <- data.frame(t(cfdata[zerocafions, ]))
cdf2 <- data.frame(t(cfdata[zerodecafions, ]))
cfdata2 <- cbind(cdf1, cfn)
cfdata3 <- cbind(cdf2, cfn)

cfmelt <- melt(cfdata2)
cfmelt2 <- melt(cfdata3)
ggplot(cfmelt, aes(x=cfn,  y=log(value))) + geom_boxplot(outlier.size=1.5,  outlier.shape=21) +
  facet_wrap(~ variable) + theme_bw()
ggplot(cfmelt2, aes(x=cfn,  y=log(value))) + geom_boxplot(outlier.size=1.5,  outlier.shape=21) +
  facet_wrap(~ variable) + theme_bw()

######################################################################################

#heatmap
#coffeemelt <- melt(allcoffees)
ggplot(cfmelt, aes(x=Var1, y=Var2, fill=log10(value))) + geom_tile() + theme_minimal() +
  scale_fill_gradient2(midpoint=4,  mid="white")
ggsave("Bar brew method vs roast.png", height= 100, width=200, units="mm")

library(latticeExtra)
#levelplot(log10(value) ~ Var1 + Var2, data=cfmelt)
