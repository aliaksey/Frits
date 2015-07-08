rm(list=ls())
# load("data/ph_raw_data.RData")
# load("alp_data.RData")
#make nice plot of  s curve
## make new plot
MAlp<-read.csv("data/CSV Results/res2_Cell_Intensity_IntegratedIntensity_ALP.csv")
for_plot<-MAlp[,c(14,12)]
##rank surfasec based on mestimation
for_plot_s<-for_plot[order(for_plot$ranksum),]
library(ggplot2)
ggplot(for_plot_s,aes(c(1:2178),mean))+geom_point(aes(size = 2, col="grey"))+
  ylab("Mean Integrated Alp Intensity per Surface")+
  xlab("Surface Rank by Ranksum")+theme_bw()+ 
  theme(legend.position = "none",text = element_text(size=20))+ 
  scale_colour_manual(values=c("grey"="grey30"))

for_plot<-MAlp[,c(13,12)]
##rank surfasec based on mestimation
for_plot_s<-for_plot[order(for_plot$ranksum),]
library(ggplot2)
ggplot(for_plot_s,aes(c(1:2178),median))+geom_point(aes(size = 2, col="grey"))+
  ylab("Median Integrated Alp Intensity per Surface")+
  xlab("Surface Rank by Ranksum")+theme_bw()+ 
  theme(legend.position = "none",text = element_text(size=20))+ 
  scale_colour_manual(values=c("grey"="grey30"))

for_plot<-MAlp[,c(14,16)]
##rank surfasec based on mestimation
for_plot_s<-for_plot[order(for_plot$mestimation),]
library(ggplot2)
ggplot(for_plot_s,aes(c(1:2178),mean))+geom_point(aes(size = 2, col="grey"))+
  ylab("Mean Integrated Alp Intensity per Surface")+
  xlab("Surface Rank by M-Estimation")+theme_bw()+ 
  theme(legend.position = "none",text = element_text(size=20))+ 
  scale_colour_manual(values=c("grey"="grey30"))

for_plot<-MAlp[,c(13,16)]
##rank surfasec based on mestimation
for_plot_s<-for_plot[order(for_plot$mestimation),]
library(ggplot2)
ggplot(for_plot_s,aes(c(1:2178),median))+geom_point(aes(size = 2, col="grey"))+
  ylab("Median Integrated Alp Intensity per Surface")+
  xlab("Surface Rank by M-Estimation")+theme_bw()+ 
  theme(legend.position = "none",text = element_text(size=20))+ 
  scale_colour_manual(values=c("grey"="grey30"))

#collapce to Image Number
###############################################################
#####################plot that used before#####################
# library(plyr)
# alp_integr_intensities<-ddply(alp_data,"ImageNumber", summarise, 
#                        AlpTrMeanIntegrativeI=mean(ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr,trimm=0.2),
#                        AlpMedianIntegrativeI=median(ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr),
#                                   FeatureIdx=unique(FeatureIdx))
# #collapce to Feature Number
# alp_integr_intensities.f<-ddply(alp_integr_intensities,"FeatureIdx", summarise, 
#      AlpTrMeanFeatureMean=mean(AlpTrMeanIntegrativeI,trimm=0.2),
#      AlpTrMeanFeatureMedian=mean(AlpMedianIntegrativeI,trimm=0.2),
#      AlpmedianFeatureMedian=median(AlpMedianIntegrativeI),
#      AlpmedianFeatureMean=mean(AlpMedianIntegrativeI))
# save(alp_integr_intensities.f,file="alp_integrated_intens.RDATA")
# hist(alp_integr_intensities.f$AlpTrMeanFeatureMedian)
# 
# library(ggplot2)
# ##another aproach without outliers
# alp_integr_intensities.f.stats<-ddply(alp_integr_intensities,"FeatureIdx", summarise, 
#                                       ymin=boxplot.stats(AlpMedianIntegrativeI)$stats[1],
#                                       lower=boxplot.stats(AlpMedianIntegrativeI)$stats[2],
#                                       middle=boxplot.stats(AlpMedianIntegrativeI)$stats[3],
#                                       upper=boxplot.stats(AlpMedianIntegrativeI)$stats[4],
#                                       ymax=boxplot.stats(AlpMedianIntegrativeI)$stats[5])
# ##ordering
# alp_integr_intensities.fs.stats<-alp_integr_intensities.f.stats[order(alp_integr_intensities.f.stats$middle),]
# alp_integr_intensities.fs.stats$FeatureIdx<-factor(alp_integr_intensities.fs.stats$FeatureIdx,
#                                                    levels=as.character(alp_integr_intensities.fs.stats$FeatureIdx))
# 
# ggplot(alp_integr_intensities.fs.stats, aes(x=FeatureIdx))+
#   geom_ribbon(data=alp_integr_intensities.fs.stats,
#               aes(group = 1,ymin=ymin, ymax=ymax,fill="blue",colour="blue"))+
#   geom_ribbon(data=alp_integr_intensities.fs.stats,
#               aes(group = 1,ymin=lower, ymax=upper,fill="yellow",colour="yellow"))+
#   geom_line(data=alp_integr_intensities.fs.stats,
#             aes(group = 1,y=middle,colour="red"))+
#   scale_fill_manual(values=c("blue"="blue","yellow"="yellow","red"="red"))+
#   scale_colour_manual(values=c("blue"="blue","yellow"="yellow","red"="red"))+
#   ylab("Median Integrated Alp Intensity per repeat")
# 
# ggplot(alp_integr_intensities.fs.stats, aes(x=FeatureIdx))+
#   geom_ribbon(data=alp_integr_intensities.fs.stats,
#               aes(group = 1,ymin=lower, ymax=upper,fill="orange",colour="orange"))+
#   geom_line(data=alp_integr_intensities.fs.stats,
#             aes(group = 1,y=middle,colour="red"))+
#   scale_fill_manual(values=c("orange"="orange","red"="red"))+
#   scale_colour_manual(values=c("orange"="orange","red"="red"))+
#   ylab("Median Integrated Alp \n Intensity per repeat")+
#   theme(axis.text.x =element_text(size  = 1,
#                                   angle = 45,
#                                   hjust = 1,
#                                   vjust = 1),legend.position="none")
# 
# ggplot(alp_integr_intensities.fs.stats, aes(FeatureIdx, 
#                                             lower=lower, upper=upper, middle=middle, ymin=ymin, ymax=ymax)) + 
#   geom_boxplot(stat="identity", fill = "#E69F00",colour = NA)+
#   geom_point(data=alp_integr_intensities.fs.stats, aes(x=c(1:2177),y=middle),
#              colour="red",shape="-",size=3)+ylab("Median Integrated Alp Intensity per repeat")+
#   theme(axis.text.x =element_text(size  = 2,
#                                   angle = 45,
#                                   hjust = 1,
#                                   vjust = 1))
# 
############################################################################
############################################################################

#save(alp_integr_intensities,alp_integr_intensities.f,file=("Alp_intensities_rank_Alex.RData"))

# plot(alp_integr_intensities.fs$AlpTrMeanFeatureMean)
# plot(alp_integr_intensities.fs$AlpTrMeanFeatureMedian)
# ##
# ##sort data based by Mhulsamn analysis
# IIALP.MH<-read.csv("data/CSV Results/res2_Cell_Intensity_IntegratedIntensity_ALP.csv")
# MH.sort<-IIALP.MH[order(IIALP.MH$mestimation),c("feature.idx","median")]
# 
# ratiorankcorr.s.mh<-merge(MH.sort,alp_integr_intensities.f,
#           by.y="FeatureIdx",by.x="feature.idx",sort=F)
# plot(ratiorankcorr.s.mh$AlpTrMeanFeatureMean)
# plot(MH.sort$median)
# 
# ggplot(ratiorankcorr.s.mh,aes(c(1:2178),AlpTrMeanFeatureMean))+
#   geom_point()+
#   geom_point()+
#   stat_smooth() 
# 
##plotting each repaeat individually
##prepare data for plots
#order data based on feature mean integrated intensities
# alp_integr_intensities.fs<-alp_integr_intensities.f[order(alp_integr_intensities.f$AlpmedianFeatureMedian),]
# #alp_integr_intensities.fs<-alp_integr_intensities.fs[c(1:100,2077:2177),]
# plot(alp_integr_intensities.fs$AlpmedianFeatureMedian)
# 
# library(reshape2)
# ratiorankcorr.s<-merge(alp_integr_intensities.fs,alp_integr_intensities,by="FeatureIdx",sort=F)
# ratiorankcorr_scurve.tr <- transform(ratiorankcorr.s[,c("FeatureIdx", "AlpMedianIntegrativeI") ], 
# FeatureIdx = factor(FeatureIdx,levels = unique(as.character(ratiorankcorr.s$FeatureIdx))))
# 
# ratiorankcorr_scurve<-melt(ratiorankcorr_scurve.tr, measure.vars ="AlpMedianIntegrativeI")
# library(ggplot2)
# ggplot(ratiorankcorr_scurve,aes(FeatureIdx,value))+
#   geom_boxplot(outlier.shape=NA,fill = "#E69F00")+
#   geom_point(data=alp_integr_intensities.fs, aes(x=c(1:2177),y=AlpmedianFeatureMedian),
#                 colour="red",shape="-")+ylab()

# ggplot(ratiorankcorr_scurve,aes(FeatureIdx,value))+
#   geom_smooth(aes(group=1),method="gam", formula=y ~ s(x, bs = "cs"),level = 0.99)
# # create plot with all statistics 
# ggplot(alp_integr_intensities.fs.stats, aes(FeatureIdx, 
#           lower=lower, upper=upper, middle=middle, ymin=ymin, ymax=ymax)) + 
#   geom_boxplot(stat="identity", fill = "#E69F00",colour = "#0072B2")+
#   geom_point(data=alp_integr_intensities.fs.stats, aes(x=c(1:2177),y=middle),
#              colour="red",shape="-")+ylab("Median Integrated Alp Intensity per repeat")
# 
# create plot with all statistics ribbon

 # #                    , 
# #                                             lower=lower, upper=upper, middle=middle, ymin=ymin, ymax=ymax)) + 
# #   geom_boxplot(stat="identity", fill = "#E69F00",colour = "#0072B2")+
# #   geom_point(data=alp_integr_intensities.fs.stats, aes(x=c(1:2177),y=middle),
# #              colour="red",shape="-")+ylab("Median Integrated Alp Intensity per repeat")
# 
# 
# 
# ggplot(alp_integr_intensities.fs.stats, aes(FeatureIdx, 
#                                             lower=lower, upper=upper, middle=middle, ymin=ymin, ymax=ymax)) + 
#   geom_boxplot(stat="identity", fill = "#E69F00",colour = "#0072B2")+
#   geom_point(data=alp_integr_intensities.fs.stats, aes(x=c(1:2177),y=middle),
#              colour="red",shape="-")+ylab("Median Integrated Alp Intensity per repeat")
# 


# create plot only with 25 75 persentile !!!!!!!!selected for paper
# 
# ##plot only 100 top and 100 bottom
# ggplot(alp_integr_intensities.fs.stats[c(1:100,2078:2177),], aes(FeatureIdx, 
#        lower=lower, upper=upper, middle=middle, ymin=ymin, ymax=ymax)) + 
#   geom_boxplot(stat="identity", fill = "#E69F00",colour = NA)+
#   geom_point(data=alp_integr_intensities.fs.stats[c(1:100,2078:2177),], 
#              aes(x=c(1:200),y=middle),
#              colour="red",shape="-")+ylab("Median Integrated Alp Intensity per repeat")
# ##plot only 100 top and 100 bottomwith all sttistics
# ggplot(alp_integr_intensities.fs.stats[c(1:100,2078:2177),], aes(FeatureIdx, 
#                                                                  lower=lower, upper=upper, middle=middle, ymin=ymin, ymax=ymax)) + 
#   geom_boxplot(stat="identity", fill = "#E69F00",colour = "#0072B2")+
#   geom_point(data=alp_integr_intensities.fs.stats[c(1:100,2078:2177),], 
#              aes(x=c(1:200),y=middle),
#              colour="red",shape="-")+ylab("Median Integrated Alp Intensity per repeat")
# 
# ##adding cours of selected hits
# surf.meta.data<-as.data.frame(cbind(FeatureIdx=hit.surf,SORT=as.character(c(1:14,"np")),Hit=
#                                        c("high","high","high","high","low","high","low","high","high","high","high","low","high","low","blank")))
# forplotcoll.fr<-merge(large.surf.meta,forplotcoll.fr.temp,sort=F)
# 
# surf.meta.data<-as.data.frame(cbind(FeatureIdx=c(1:2177),Hit=as.character(rep("nonsel",2177))),stringsAsFactors=F)
# 
# surf.meta.data[surf.meta.data$FeatureIdx%in%c(11,42,46,1050,2114,1018,1673,2150,443,689),"Hit"]<-"High"
# surf.meta.data[surf.meta.data$FeatureIdx==2177,"Hit"]<-"NP"
# surf.meta.data[surf.meta.data$FeatureIdx%in%c(1147,229,864,494),"Hit"]<-"Low"
# 
# alp_integr_intensities.fs.stats.m<-merge(alp_integr_intensities.fs.stats,surf.meta.data,by="FeatureIdx",sort=F)
# 
# ggplot(alp_integr_intensities.fs.stats.m, aes(FeatureIdx, 
#                                               lower=lower, upper=upper, middle=middle, 
#                                               ymin=ymin, ymax=ymax)) + 
#   geom_boxplot(stat="identity", colour = NA,aes(fill = factor(Hit)))+
#   geom_point(data=alp_integr_intensities.fs.stats, aes(x=c(1:2177),y=middle),
#              colour="red",shape="-")+ylab("Median Integrated Alp Intensity per repeat")
# 
# # 
# # stat_sum_df <- function(fun, geom="crossbar", ...) { 
# #   stat_summary(fun.data=fun, colour="blue", geom=geom, width=0.2, ...) 
# # } 
# # 
# # ggplot(alp_integr_intensities.fs,aes(c(1:2177),AlpTrMeanFeatureMedian))+
# #   geom_point()+
# #   stat_sum_df("mean_cl_normal", geom = "smooth") 
# # 
# # ggplot(alp_integr_intensities.fs,aes(c(1:2177),AlpTrMeanFeatureMedian))+
# #   geom_point()+
# #   geom_smooth()
# # 
# # ggplot(alp_integr_intensities.fs,aes(c(1:2177),AlpTrMeanFeatureMedian))+
# #   geom_point()+
# #   stat_summary(fun.y = 'mean', colour = 'blue', geom = 'line')
# # stat_summary(fun.data = 'mean_sdl', geom = 'ribbon', alpha = 0.2)
# # 
# # ggplot(alp_integr_intensities.fs,aes(c(1:2177),AlpTrMeanFeatureMedian))+
# #   geom_point()+
# #   geom_point()+
# #   stat_smooth() 
# ##calculate kruskal wallis statistics for different intensities
# 
