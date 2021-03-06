rm(list=ls())
load("data//ph_raw_data.RData")
alp_data<-read.csv("data/data_all_alp_data_per_object_old.csv",stringsAsFactors=F,header = F)
alp_data.heads<-read.csv("data/data_all_alp_data_per_object_old_heads.csv",stringsAsFactors=F)
if(all(alp_data.heads==alp_data[1:1000,])) colnames(alp_data)<-colnames(alp_data.heads)
colnames(alp_data)
alp_data<-merge(alp_data,image.data[,c("ImageNumber", "FeatureIdx")], by="ImageNumber")
save(alp_data,file="alp_data.RData")
