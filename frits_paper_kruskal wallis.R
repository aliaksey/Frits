rm(list=ls())
load("data/ph_raw_data.RData")
load("alp_data.RData")
#make nice plot of  s curve
#collapce to Image Number
alp_data$ALPCytoplasm_Intensity_RelativeIntensity_ALP4Corr<-
  alp_data$ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr/
  alp_data$ALPCytoplasm_Intensity_IntegratedIntensity_Actin4Corr
library(plyr)
alp_intensities<-ddply(alp_data,"ImageNumber", summarise, 
  AlpMedianIntegrativeI=median(ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr),
  AlpMedianMaxI=median(ALPCytoplasm_Intensity_MaxIntensity_ALP4Corr),
  AlpMedianRelI=median(ALPCytoplasm_Intensity_RelativeIntensity_ALP4Corr),
  FeatureIdx=unique(FeatureIdx))
                              
#   FeatureIdx=unique(FeatureIdx))
# #collapce to Feature Number
# alp_integr_intensities.f<-ddply(alp_integr_intensities,"FeatureIdx", summarise, 
#                                 AlpmedianFeatureMean=mean(AlpMedianIntegrativeI))
##calculating Kruskal wallis statistics

pvintg<-kruskal.test(AlpMedianIntegrativeI ~ FeatureIdx, data = alp_intensities)$p.value
pvmax<-kruskal.test(AlpMedianMaxI ~ FeatureIdx, data = alp_intensities)$p.value
pvrel<-kruskal.test(AlpMedianRelI ~ FeatureIdx, data = alp_intensities)$p.value

kwplot<-as.data.frame(cbind(Intensity=c("Integrative","Maximum","Relative"),
                            Pvalue=c(as.numeric(pvintg),as.numeric(pvmax),
                                     as.numeric(pvrel))),stringsAsFactors=F)
kwplot$Pvaluelog<-abs(log10(as.numeric(kwplot$Pvalue)))

qplot(factor(kwplot$Intensity), kwplot$Pvaluelog, geom="bar", stat="identity")+ coord_flip()+
  geom_hline(yintercept = abs(log10(0.5)), colour="red")+ylab("-log10 p-value Kruskal wallis test")+
  xlab("ALP Intensity")

##krusakl wallis for all dat
pvintg.m<-kruskal.test(ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr ~ FeatureIdx, data = alp_data)$p.value
pvmax.m<-kruskal.test(ALPCytoplasm_Intensity_MaxIntensity_ALP4Corr ~ FeatureIdx, data = alp_data)$p.value
pvrel.m<-kruskal.test(ALPCytoplasm_Intensity_RelativeIntensity_ALP4Corr ~ FeatureIdx, data = alp_data)$p.value

kwplot<-as.data.frame(cbind(Intensity=c("Integrative","Maximum","Relative"),
                            Pvalue=c(as.numeric(pvintg.m),as.numeric(pvmax.m),
                                     as.numeric(pvrel.m))),stringsAsFactors=F)
kwplot$Pvaluelog<-abs(log10(as.numeric(kwplot$Pvalue)))

qplot(factor(kwplot$Intensity), kwplot$Pvaluelog, geom="bar", stat="identity")+ coord_flip()+
  geom_hline(yintercept = abs(log10(0.5)), colour="red")+ylab("-log10 p-value Kruskal wallis test")+
  xlab("ALP Intensity")