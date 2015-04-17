rm(list=ls())
load("Alp_intensities_rank_Alex.RData")
integrated_intensity.temp<-read.csv("data/CSV Results/res2_Cell_Intensity_IntegratedIntensity_ALP.csv",stringsAsFactors=F)
max_intensity.temp<-read.csv("data/CSV Results/res2_Cell_Intensity_MaxIntensity_ALP.csv",stringsAsFactors=F)
relative_intensity.temp<-read.csv("data/CSV Results/res2_Cell_Intensity_RelativeIntegratedIntensity_ALP.csv",stringsAsFactors=F)
##remove nonrelevant data
integrated_intensity<-integrated_intensity.temp[,-c(3:10,17:57)]
max_intensity<-max_intensity.temp[,-c(2:11,17:57)]
relative_intensity<-relative_intensity.temp[,-c(2:11,17:57)]
##fix names
colnames(integrated_intensity)[!colnames(integrated_intensity)%in%c("internal.idx",
                "feature.idx","number.of.outliers")]<-paste("Intg_int",colnames(integrated_intensity)[!
              colnames(integrated_intensity)%in%c("internal.idx","feature.idx",
                                                  "number.of.outliers")],sep="_")

colnames(max_intensity)[!colnames(max_intensity)%in%c("internal.idx",
          "feature.idx")]<-paste("Max_int",colnames(max_intensity)[!
   colnames(max_intensity)%in%c("internal.idx","feature.idx")],sep="_")

colnames(relative_intensity)[!colnames(relative_intensity)%in%c("internal.idx",
 "feature.idx")]<-paste("Rel_int",colnames(relative_intensity)[!
    colnames(relative_intensity)%in%c("internal.idx","feature.idx")],sep="_")
##join evrything in one table
alp_mark_rank<-merge(integrated_intensity, merge(max_intensity,relative_intensity))

#library(GGally)
# ggpairs(integrated_intensity,columns=3:7)
# ggpairs(max_intensity,columns=3:7)
# ggpairs(relative_intensity,columns=3:7)

##show expression of alp on hit surfaces

hit.surf<-c(11,42,46,1050,1147,2114,229,1018,1673,2150,443,864,689,494,2177)
forplotcoll.fr.temp<-alp_integr_intensities[alp_integr_intensities$FeatureIdx%in%hit.surf,]
#alp_integr_intensities.fs<-alp_integr_intensities.f[order(-alp_integr_intensities.f$AlpTrMeanFeatureMean),]
#forplotcoll.fr<-merge(alp_integr_intensities.fs,forplotcoll.fr.temp, sort=F, all=F)

large.surf.meta<-as.data.frame(cbind(FeatureIdx=hit.surf,SORT=as.character(c(1:14,"np")),Hit=
  c("high","high","high","high","low","high","low","high","high","high","high","low","high","low","blank")))
forplotcoll.fr<-merge(large.surf.meta,forplotcoll.fr.temp,sort=F)
###Trimmed mean
forplottransf.fr <- transform(forplotcoll.fr[,c("SORT", "AlpTrMeanIntegrativeI","Hit")], SurfaceNumber = factor(SORT, 
                levels = unique(as.character(forplotcoll.fr$SORT))))
library(reshape2)
forplottransfmelt.fr<-melt(forplottransf.fr, measure.vars ="AlpTrMeanIntegrativeI")
library(ggplot2)
ggplot(forplottransfmelt.fr, aes(SORT, value, fill=Hit))+
  geom_boxplot(outlier.shape=NA)+theme(legend.position="none")+geom_jitter(size=3)+
  ylab("Robust Mean Integrated Alp Intensity (a.u.)")+xlab("Surface Number")+ 
  scale_x_discrete(limits=large.surf.meta$SORT)

###Median
forplottransf.fr <- transform(forplotcoll.fr[,c("SORT", "AlpMedianIntegrativeI","Hit")], SurfaceNumber = factor(SORT, 
                                                                                                                levels = unique(as.character(forplotcoll.fr$SORT))))
library(reshape2)
forplottransfmelt.fr<-melt(forplottransf.fr, measure.vars ="AlpMedianIntegrativeI")
library(ggplot2)
ggplot(forplottransfmelt.fr, aes(SORT, value, fill=Hit))+
  geom_boxplot(outlier.shape=NA)+theme(legend.position="none")+geom_jitter(size=3)+
  ylab("Median Integrated Alp Intensity (a.u.)")+xlab("Surface Number")+ 
  scale_x_discrete(limits=large.surf.meta$SORT)
#####Some extra plots 
# library(doBy)
# icam_freq_stat <- summaryBy(Ratio ~ FeatureIdx, forplotcoll.fr,order=F, FUN = c(mean, sd))
# 
# SD.quartile <- function(x){
#   mean.temp<-mean(x)
#   sd1.temp<-sd(x)
#   n.temp<-length(x)
#   CI.int<-0.995
#   error.temp <- qnorm(CI.int)*sd1.temp/sqrt(n.temp)
#   out <- c(mean.temp-error.temp,mean.temp,mean.temp+error.temp)
#   names(out) <- c("ymin","y","ymax")
#   return(out) 
# }
# library(ggplot2)
# #plot dots + confidence interval
# ggplot(forplottransfmelt.fr, aes(FeatureIdx, value,group=FeatureIdx,colour=FeatureIdx))+
#   stat_summary(fun.data = SD.quartile, geom = "pointrange",color = "Blue")+
#   geom_jitter( size = 4)  
# 
# 
# ggplot(forplottransfmelt.fr, aes(FeatureIdx, y = value, fill=FeatureIdx))+
#   geom_boxplot(outlier.shape=NA)+theme(legend.position="none")+geom_jitter(size=3)+
#   stat_summary(fun.data = SD.quartile, geom = "errorbar",color = "red")+geom_smooth(aes(group = 1),method="glm",colour="grey50")+
#   ylab("Robust Mean Integrated Alp Intensity (a.u.)")+xlab("Feature Index")
# 
# ggplot(forplottransfmelt.fr, aes(FeatureIdx, y = value, fill=FeatureIdx))+
#   geom_boxplot(outlier.shape=NA)+theme(legend.position="none")+geom_jitter(size=3)+
#   stat_summary(fun.data = SD.quartile, geom = "errorbar",color = "red")+geom_smooth(aes(group = 1),method="glm",colour="grey50")+
#   ylab("Robust Mean Integrated Alp Intensity (a.u.)")+xlab("Feature Index")+theme_bw()
# 
# 
# ##making fold change plot
# forplotcoll.fr.temp$IntenFCh<-forplotcoll.fr.temp$AlpTrMeanIntegrativeI/
#   mean(forplotcoll.fr.temp[forplotcoll.fr.temp$FeatureIdx==2177,"AlpTrMeanIntegrativeI"])
# 
# forplotcoll.fr<-merge(alp_integr_intensities.fs,forplotcoll.fr.temp, by="FeatureIdx",sort=F, all=F)
# forplottransf.fr <- transform(forplotcoll.fr[,c("FeatureIdx", "IntenFCh")], FeatureIdx = factor(FeatureIdx, 
#          levels = unique(as.character(forplotcoll.fr$FeatureIdx))))
# library(reshape2)
# forplottransfmelt.fr<-melt(forplottransf.fr, measure.vars ="IntenFCh")
# #plot dots + SD
# 
# SD.quartile2 <- function(x){
#   mean.temp<-mean(x)
#   sd1.temp<-sd(x)
#   n.temp<-length(x)
#   CI.int<-0.995
#   error.temp <- qnorm(CI.int)*sd1.temp/sqrt(n.temp)
#   out <- c(mean.temp-sd1.temp,mean.temp,mean.temp+sd1.temp)
#   names(out) <- c("ymin","y","ymax")
#   return(out) 
# }
# 
# dodge <- position_dodge(width=0.9)
# ggplot(forplottransfmelt.fr, aes(FeatureIdx, value,group=FeatureIdx,fill=FeatureIdx))+
#   stat_summary(fun.y="mean", geom="bar")+ 
#   stat_summary(fun.data = SD.quartile2, geom = "errorbar",color = "Blue")+ylab("Fold Chnge Mean Integrated Intensity")
# ggplot(forplottransfmelt.fr, aes(FeatureIdx, value,group=FeatureIdx,fill=FeatureIdx))+
#   stat_summary(fun.y="mean", geom="bar")+ 
#   stat_summary(fun.data = SD.quartile, geom = "errorbar",color = "red")+ylab("Fold Chnge Mean Integrated Intensity")

######Compare ranking for ALP max Integrated and Normolised
# ##plot scurves
# ##function to make scurve plot
# scurveplot<-function(sortvar,plotvar,data.tb){
# ##sort data based on variable
# data_sorted_intg<-data.tb[order(data.tb[,sortvar]),
#              c("feature.idx",sortvar,plotvar)]
# #prepare data for plots
# library(reshape2)
# ratiorankcorr_scurve.tr <- transform(data_sorted_intg, 
# FeatureIdx = factor(feature.idx,levels = unique(as.character(data_sorted_intg$feature.idx))))
# 
# ratiorankcorr_scurve<-melt(ratiorankcorr_scurve.tr, measure.vars =plotvar)
# library(ggplot2)
# gplot<-ggplot(ratiorankcorr_scurve,aes(FeatureIdx,value))+
#   geom_point()+geom_smooth(aes(group=1), method="glm")+ylab(plotvar)
# print(gplot)
# }
# 
# ##making plots
# #based on median
# scurveplot("Intg_int_median","Intg_int_median",alp_mark_rank)
# scurveplot("Intg_int_mestimation","Intg_int_median",alp_mark_rank)
# scurveplot("Max_int_median","Max_int_median",alp_mark_rank)
# scurveplot("Max_int_mestimation","Max_int_median",alp_mark_rank)
# scurveplot("Rel_int_median","Rel_int_median",alp_mark_rank)
# scurveplot("Rel_int_mestimation","Rel_int_median",alp_mark_rank)
# 
# #  plotting u test values
# 
# summary(alp_mark_rank[,c("Intg_int_ranksum","Max_int_ranksum" ,"Rel_int_ranksum" )])
# alp_bar_plot<-melt(alp_mark_rank, measure.vars =c("Intg_int_ranksum",
#                  "Max_int_ranksum" ,"Rel_int_ranksum" ))
# library(ggplot2)
# ggplot(alp_bar_plot,aes(factor(variable),value))+
# geom_bar(stat="identity")  
