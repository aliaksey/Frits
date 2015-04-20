featind.temp<-read.csv("data/indexes_and_features.csv")
featind<-featind.temp[,-c(1,3:13)]
alp_model.temp<-read.csv("data/CSV Results/res2_Cell_Intensity_IntegratedIntensity_ALP.csv",stringsAsFactors=F)
featind.2<-alp_model.temp[,-c(1,3:19)]
##make some checks
all(colnames(featind)==colnames(featind.2))
all(featind[order(featind$feature.idx),]==featind.2[order(featind.2$feature.idx),])

##perform some analysis
summary(featind$TriArea)

###compute new variables
colnames(featind)

##sum of absolute primitive areas
featind$SumArea<-featind$TriArea+featind$CircArea+featind$LineArea
summary(featind$SumArea)

##number of patterns per surface
featind$PatNum<-280^2/featind$FeatSize^2
summary(featind$PatNum)

##density of patterns on surface
featind$PatDen<-(featind$PatNum*featind$SumArea)/280^2
summary(featind$PatDen)
featind$PatDen2<-(featind$PatNum*(featind$FCP*(featind$FeatSize^2)))/280^2
summary(featind$PatDen2)

##actual feature area
featind$PatAr<-featind$FCP*(featind$FeatSize^2)
summary(featind$PatAr)
plot(featind$PatAr,featind$FeatSize^2)
##actual empty space area
featind$SpcAr<-(1-featind$FCP)*(featind$FeatSize^2)
summary(featind$SpcAr)
plot(featind$SpcAr,featind$PatAr)

##Difference between pattern area and sum primitives areas
featind$PatArDif<-featind$PatAr-featind$SumArea
summary(featind$PatArDif)

##total Number of Primitives
featind$TNP<-featind$NumCirc+featind$NumTri+featind$NumLine
summary(featind$TNP)

##Find features that consists fromdiscrete patterns
featind$PatAr2d<-NA
featind[featind$PatArDif>0,"PatAr2d"]<-featind[featind$PatArDif>0,"PatAr"]/featind[featind$PatArDif>0,"TNP"]
featind[featind$PatArDif<=0,"PatAr2d"]<-featind[featind$PatArDif<=0,"PatAr"]

featind$PatAr22d<-NA
featind[featind$PatArDif>0,"PatAr22d"]<-featind[featind$PatArDif>0,"PatAr"]/(featind[featind$PatArDif>0,"TNP"]/2)
featind[featind$PatArDif<=0,"PatAr22d"]<-featind[featind$PatArDif<=0,"PatAr"]

summary(featind$PatAr2d)
summary(featind$PatAr22d)

##spacing
featind$Spc<-sqrt(featind$FeatSize^2-featind$SumArea)
summary(featind$Spc)

featind$Spc2<-sqrt(featind$FeatSize^2-(featind$FeatSize^2*featind$FCP))
summary(featind$Spc2)

featind$Spc3<-sqrt((featind$FeatSize^2-(featind$FeatSize^2*featind$FCP))/2)
summary(featind$Spc3)

##Number of features per 

##replacing blank
featind[featind$feature.idx==2177,]<-0
any(is.na(featind))
##

##make corrplot
library(corrplot)
corrplot(cor(featind[,-1]))


plot(featind$PatAr2d,featind$FeatSize)
plot(featind$PatAr22d,featind$FeatSize)
plot(featind$PatArDif,featind$FeatSize)
plot(featind$PatAr,featind$FeatSize)
plot(featind$PatAr,featind$SumArea)
cor(featind$PatAr,featind$SumArea)

hist(featind$FCP)


plot(featind$Spc,featind$FeatSize)
plot(featind$Spc2,featind$FeatSize)
plot(featind$Spc,featind$Spc3)

save(featind,file="SurfaceFeaturesPlus.RData")
