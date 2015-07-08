rm(list=ls())
library(caret)
library(corrplot)
library(randomForest)
library(rpart)
library(partykit)
library(gplots)
set.seed(683475823)
alp_model.temp<-read.csv("data/CSV Results/res2_Cell_Intensity_IntegratedIntensity_ALP.csv",stringsAsFactors=F)
##cleaning variables
alp_model<-alp_model.temp[,-c(1,3:11,17:19)]

##loading all per cell data
#load("D:/projects/phenome_project/phenome/Cell all data & ground truth scaled.RData")
load("data/ph_raw_data.RData")

#attach alp rank to morphological data
alp_model_morph.temp<-merge(alp_model,cell.ftrs,by.x="feature.idx",by.y="FeatureIdx")
#Selecting features for anlaysis
alp_model_morph<-alp_model_morph.temp[,c("feature.idx","median","ImageNumber","ObjectNumber",
                                         colnames(alp_model_morph.temp)[grepl("Cells_AreaShape_",colnames(alp_model_morph.temp))&
                                                                          !grepl("_Center_",colnames(alp_model_morph.temp))])]
colnames(alp_model_morph)

###get rid of highly correlated features
# corMatips <- cor(alp_model_morph[,-c(1:4)], method="spearman")
# corrplot(corMatips)

# hCorr <- findCorrelation(corMatips, 0.75)
# alp_model_f<-alp_model_morph[,!colnames(alp_model_morph)%in%colnames(corMatips[,hCorr])]

alp_model_f<-alp_model_morph

colnames(alp_model_f)

corMatips.f <- cor(alp_model_f[,-c(1:4)],method="spearman")
corrplot(corMatips.f)
alp_model_f<-alp_model_morph

##rename rows
# alp_model_f[alp_model_f$feature.idx==2177,"feature.idx"]<-paste(2177,c(1:2),sep="_")
alp_model_f<-alp_model_f[alp_model_f$feature.idx!=2177,]
rownames(alp_model_f)<-paste(alp_model_f$feature.idx,rownames(alp_model_f),sep="_")
# #add new ranking
# load("Alp_intensities_rank_Alex.RData")
# alp_model_f<-merge(alp_model_f,alp_integr_intensities.f,by.x="feature.idx",by.y="FeatureIdx")
# alp_model_f$median.new<-alp_model_f$AlpTrMeanFeatureMean
# alp_model_f<-alp_model_f[,!colnames(alp_model_f)%in%c("AlpTrMeanFeatureMean","AlpTrMeanFeatureMedian")]
#Selecting hit 100 from both sides
alp_model_f.s<-alp_model_f[order(alp_model_f$median),]
plot(alp_model_f.s$median)
##select 200 features
top.feat<-tail(unique(alp_model_f.s$feature.idx),100)
bottom.feat<-head(unique(alp_model_f.s$feature.idx),100)

alp_model_f.ss.t<-alp_model_f.s[alp_model_f.s$feature.idx%in%top.feat,]
alp_model_f.ss.t$Class<-"pos"
alp_model_f.ss.b<-alp_model_f.s[alp_model_f.s$feature.idx%in%bottom.feat,]
alp_model_f.ss.b$Class<-"neg"
alp_model_f.ss<-rbind(alp_model_f.ss.b,alp_model_f.ss.t)
plot(alp_model_f.ss.t$median)
plot(alp_model_f.ss.b$median)
plot(alp_model_f.ss$median)
##find most severe outliers (boxplot rule)
#alp_model_f.sso<-c()
#library(plyr)
rm(list=c("rslt","rsltbb","rslta"))
for(i in unique(alp_model_f.ss[,"feature.idx"])){
  temp<-alp_model_f.ss[alp_model_f.ss$feature.idx==i,]
  ###filter based on all features cell number
  for(f in c(5:(length(temp)-1))){
    temp2<-temp[,c(3,4,f)]
    #temp2$CellIdx<-paste(temp2$ImageNumber,temp2$ObjectNumber,sep="_")
    lbnda<-as.numeric(quantile(temp2[,3], probs = 0.25))
    ubnda<-as.numeric(quantile(temp2[,3], probs = 0.75))
    iuda<-ubnda-lbnda
    rslta<-temp2[temp2[,3]<(ubnda+1.5*iuda)&
                   temp2[,3]>(lbnda-1.5*iuda),]
    if (length(rslta[,1])==0) break
    if(!exists("rsltbb"))rsltbb<-rslta else rsltbb<-merge(rslta, rsltbb, 
                                                          by=c("ImageNumber","ObjectNumber"))
    
  }
  if(!exists("rslt"))rslt<-rsltbb else rslt<-rbind(rsltbb, rslt)
  if(f==(length(temp)-1)) rm("rsltbb")
}
alp_model_f.sso<-merge(alp_model_f.ss[,c(1:3,length(alp_model_f.ss))],rslt,by="ImageNumber") 
alp_model_f.sso$Class<-as.factor(alp_model_f.sso$Class)
##take median per repeat
library(plyr)
alp_model_f.sso2<-ddply(alp_model_f.sso,"ImageNumber", numcolwise(median))
# alp_model_f.sso3<-merge(alp_model_f.sso2,alp_model_f.sso[,c(1,4)],by="ImageNumber")
#median per feature
alp_model_f.sso2.f<-ddply(alp_model_f.sso2,"feature.idx", numcolwise(median))

alp_model_f.sso3.f<-merge(alp_model_f.sso2.f,unique(alp_model_f.sso[,c(2,4)]),
                          by="feature.idx")

alp_model_f.sss<-alp_model_f.sso3.f[,-c(1:4)]
##creating model

##selecting samples for training and testing
data_for_model<-alp_model_f.sss
colnames(data_for_model)<-gsub("Cells_AreaShape_", "",colnames(data_for_model))
colnames(data_for_model)
#data_for_model<-alp_model_f.sss[,colnames(alp_model_f.sss2)!="median"]
##partition basd on feature
feature_class<-unique(alp_model_f.sso3.f[,c(1,length(alp_model_f.ss))])
class_data_feat<-feature_class[,"Class"]
FeatureInTrain<-createDataPartition(class_data_feat, p =3/4, list = FALSE)
inTrain<-alp_model_f.sso3.f$feature.idx%in%feature_class[FeatureInTrain,"feature.idx"]
forTraining <- data_for_model[inTrain,]
#forTraining <- ipshitsmodel
forTrainingX <- forTraining[, names(forTraining) != "Class"]
#reate testing set for features selection
forTesting <- data_for_model[!(inTrain),]
#forTesting <- ipshitsmodel

#############rpart analysis##########
#     rpart_training <- rpart(Class~.,  method="class", data=data_for_model)
#     #plot(rpart_training)
#     #text(rpart_training)
#     ##polt as party object
#     rpart1a <- as.party(rpart_training)
#     plot(rpart1a, main="Pruned CART classyfication tree for OCT4 hits")
#     
##rpart on unchanged data

rpart_training2 <- rpart(Class~.,  method="class", data=data_for_model)
#plot(rpart_training)
#text(rpart_training)
##polt as party object
rpart1a2 <- as.party(rpart_training2)
plot(rpart1a2)#, main="Pruned CART classyfication tree for OCT4 hits")


# rpartPred <- predict(rpart_training, forTesting, type = "class")
# confusionMatrix(rpartPred, forTesting$Class)
#______________________________________________________________________________________ 
##tunning the model
cvCtrl <- trainControl(method = "repeatedcv", repeats = 10,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,savePred=T,returnResamp="final")
# #rpart
# rpartTune <- train(Class ~ ., data = forTraining, method = "rpart",
#                    tuneLength = 10,
#                    metric = "ROC",
#                    trControl = cvCtrl)
# plot(rpartTune)
# predictors(rpartTune)
# plot(varImp(rpartTune),Positive=5,cex=4,pch=16,
#      main="Feature importance for CART method")
# varImp_re<-varImp(rpartTune)
# row.names(varImp_re$importance)[varImp_re$importance>0]
# plot.train(rpartTune)
# #print.train(rpartTune)
# plot(rpartTune, scales = list(x = list(log = 10)))
# 
# rpartPred2 <- predict(rpartTune, forTesting)
# confusionMatrix(rpartPred2, forTesting$Class)
# table(forTesting$Class)
# 
# rpartProbs <- predict(rpartTune, forTesting, type = "prob")
# head(rpartProbs)
# library(pROC)
# rpartROC <- roc(predictor = rpartProbs$Positive,
#                 response = forTesting$Class,
#                 levels = rev(levels(forTesting$Class)))
# 
# plot(rpartROC, type = "S", print.thres = .5)
# rpartROC
# ##svm
# svmTune <- train(Class ~ ., data = forTraining,
#                  method = "svmRadial",
#                  tuneLength = 10,
#                  preProc = c("center", "scale"),
#                  metric = "ROC",
#                  trControl = cvCtrl)#,
# #classProbs =  TRUE)
# #svmTune
# 
# #predictors(svmTune)
# plot(varImp(svmTune),
#      main="Feature importance for SVM method")
# svmTune$finalModel
# plot(svmTune, metric = "ROC", scales = list(x = list(log =2)))
# svmPred <- predict(svmTune, forTesting[, names(forTesting) != "Class"])
# confusionMatrix(svmPred, forTesting$Class)
# 
# svmProbs <- predict(svmTune, forTesting[, names(forTesting) != "Class"],
#                     type = "prob")
# head(svmProbs)
# svmROC <- roc(predictor = svmProbs$Positive,
#               response = forTesting$Class,
#               levels = rev(levels(forTesting$Class)))
# plot(svmROC, type = "S", print.thres = .5)
# #svmROC
# 
# # ##logit
# logitTune <- train(x = forTrainingX,
#                    y = forTraining$Class,
#                    method = "glm",
#                    tuneLength = 10,
#                    family = binomial(link = "logit"),
#                    metric = "ROC",
#                    trControl = cvCtrl)
# summary(logitTune)
# 
# #predictors(logitTune)
# plot(varImp(logitTune),Positive=17,cex=4,pch=16,cex.axis=30,
#      main="Feature importance for Logit method")
# logitTune$finalModel
# logitPred <- predict(logitTune, forTesting[, names(forTesting) != "Class"])
# confusionMatrix(logitPred, forTesting$Class)
# 
# logitProbs <- predict(logitTune, forTesting[, names(forTesting) != "Class"],
#                       type = "prob")
# head(logitProbs)
# logitROC <- roc(predictor = logitProbs$Positive,
#                 response = forTesting$Class,
#                 levels = rev(levels(forTesting$Class)))
# plot(logitROC, type = "S", print.thres = .5)
# 
# ##nb
# nbTune <- train(x = forTrainingX,
#                 y = forTraining$Class,
#                 method = "nb",
#                 tuneLength = 10,
#                 metric = "ROC",
#                 trControl = cvCtrl)
# summary(nbTune)
# 
# predictors(nbTune)
# plot(varImp(nbTune),
#      main="Feature importance for NB method")
# nbTune$finalModel
# nbPred <- predict(nbTune, forTesting[, names(forTesting) != "Class"])
# confusionMatrix(nbPred, forTesting$Class)
# 
# nbProbs <- predict(nbTune, forTesting[, names(forTesting) != "Class"],
#                    type = "prob")
# head(nbProbs)
# nbROC <- roc(predictor = nbProbs$Positive,
#              response = forTesting$Class,
#              levels = rev(levels(forTesting$Class)))
# plot(nbROC, type = "S", print.thres = .5)

##rf
library(doParallel)
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
rfTune <- train(x = forTrainingX,
                y = forTraining$Class,
                method = "rf",
                tuneLength = 10,
                allowParallel=TRUE,
                metric = "ROC",
                trControl = cvCtrl)
stopCluster(cl)
registerDoSEQ()


summary(rfTune)

predictors(rfTune)
plot(varImp(rfTune), top=20,cex=1.5,cex.lab=3)

rfTune$finalModel
rfPred <- predict(rfTune, forTesting[, names(forTesting) != "Class"])
confusionMatrix(rfPred, forTesting$Class)

rfProbs <- predict(rfTune, forTesting[, names(forTesting) != "Class"],
                   type = "prob")
head(rfProbs)
rfROC <- roc(predictor = rfProbs$pos,
             response = forTesting$Class,
             levels = rev(levels(forTesting$Class)))
plot(rfROC, type = "S", print.thres = 0.356,col = "red",lwd=2.5,
     ylab="True Positive Rate(Sensitivity)",
     xlab="False Positive rate (1-Specificity)")
legend(x=0.6,y=0.2,
       legend=c(paste(" AUC=","0.77"),
                paste("Accuracy=","0.7" )),box.lwd = 0,box.col = "white",bg = "white")

library(ROCR)
predrf <- prediction(rfProbs$pos, forTesting$Class)
perfrf <- performance( predrf, "tpr", "fpr")
plot( perfrf, col="red",lty=1, lwd=3)
abline(a=0,b=1,col="grey",lwd=1,lty=3)  
auc<-performance(predrf,"auc")
auc <- unlist(slot(auc, "y.values"))
auct <- paste(c("AUC  = "),auc,sep="")
legend(x=0.6,y=0.3,
       legend=c(paste(" AUC=","0.78"),
                paste("Accuracy=","0.72" )),box.lwd = 0,box.col = "white",bg = "white")

auc
# 
# ##LDA
# ldaTune <- train(x = forTrainingX,
#                  y = forTraining$Class,
#                  method = "lda",
#                  tuneLength = 10,
#                  metric = "ROC",
#                  trControl = cvCtrl)
# summary(ldaTune)
# 
# predictors(ldaTune)
# plot(varImp(ldaTune),Positive=17,cex=4,pch=16,cex.axis=30,
#      main="Feature importance for lda method")
# ldaTune$finalModel
# ldaPred <- predict(ldaTune, forTesting[, names(forTesting) != "Class"])
# confusionMatrix(ldaPred, forTesting$Class)
# 
# ldaProbs <- predict(ldaTune, forTesting[, names(forTesting) != "Class"],
#                     type = "prob")
# head(ldaProbs)
# ldaROC <- roc(predictor = ldaProbs$Positive,
#               response = forTesting$Class,
#               levels = rev(levels(forTesting$Class)))
# plot(ldaROC, type = "S", print.thres = .5)
# 
# ##QDA
# qdaTune <- train(x = forTrainingX,
#                  y = forTraining$Class,
#                  method = "qda",
#                  tuneLength = 10,
#                  metric = "ROC",
#                  trControl = cvCtrl)
# summary(qdaTune)
# 
# predictors(qdaTune)
# plot(varImp(qdaTune),Positive=17,cex=4,pch=16,cex.axis=30,
#      main="Feature importance for qda method")
# qdaTune$finalModel
# qdaPred <- predict(qdaTune, forTesting[, names(forTesting) != "Class"])
# confusionMatrix(qdaPred, forTesting$Class)
# 
# qdaProbs <- predict(qdaTune, forTesting[, names(forTesting) != "Class"],
#                     type = "prob")
# head(qdaProbs)
# qdaROC <- roc(predictor = qdaProbs$Positive,
#               response = forTesting$Class,
#               levels = rev(levels(forTesting$Class)))
# plot(qdaROC, type = "S", print.thres = .5)
# 
# 
# ##compare 
# plot(rpartROC, type = "S",col = "red", main="comparison 3 fitted models")#,print.thres = .5
# plot(svmROC, add = TRUE, col = "blue")#,print.thres = .5
# plot(logitROC, add = TRUE, col = "yellow")
# plot(nbROC, add = TRUE, col = "grey")
# plot(rfROC, add = TRUE, col = "green")
# plot(ldaROC, add = TRUE, col = "black")
# plot(qdaROC, add = TRUE, col = "purple")
# legend(x=0.45,y=0.5, legend=c(paste("CART: AUC=",as.numeric(rpartROC$auc)),
#                               paste("SVM: AUC=",as.numeric(svmROC$auc)),
#                               paste("Logit: AUC=",as.numeric(logitROC$auc)),
#                               paste("NB: AUC=",as.numeric(nbROC$auc)),
#                               paste("RF: AUC=",as.numeric(rfROC$auc)),
#                               paste("LDA: AUC=",as.numeric(ldaROC$auc)),
#                               paste("QDA: AUC=",as.numeric(qdaROC$auc))),
#        lty=1, col=c("red","blue","yellow","grey","green","black","purple"), bty='o', cex=1)
# 
# save(list=ls(),file="classyfication cell per feature.RDATA")
# ###find most important parameters based on methods an
#For SVM classification models, the default behavior 
#is to compute the area under the ROC curve. 
# plot(varImp(svmTune))
# plot(varImp(svmTune,useModel = T,scale = FALSE))
# #area under the ROC curve for each predictor
# rocimpsvm<-filterVarImp(x = forTraining[, -ncol(forTraining)], y = forTraining$Class)
# rocimpsvm[order(-rocimpsvm[,2]),]
# 
# ## plot alp intensity values base on that parameters
 library(ggplot2)
ggplot(alp_model_f.sso3.f,aes(Cells_AreaShape_Solidity,
                              Cells_AreaShape_Compactness,shape=Class,colour=median))+geom_point()
ggplot(data_for_model,aes(Compactness,
                          Solidity,colour=as.factor(Class)))+geom_point()+
  theme_bw()+
  scale_color_discrete(name ="Surface Class for \n Alp expression", 
                       labels=c("Negative","Positive"))+ theme(legend.position = c(.8, .8))
# 
# create per cell plot
# predare data select image that passed outlier removal
# select only relevant features
load("alp_data.RData")
# library(ggplot2)
colnames(alp_data)<-gsub("Cells_AreaShape_", "",colnames(alp_data))

cell_plot<-alp_data[,c("ImageNumber",
                       "ObjectNumber","Compactness",
                       "Solidity",
                       "ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr")]
##select images and objects passed QC
cell_plot.f<-unique(merge(cell_plot,
                          alp_model_f.sso[,c(1:5)],by=c("ImageNumber","ObjectNumber")))
cell_plot.f$ALP_Int_scale<-scale(cell_plot.f$ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr)

# ggplot(cell_plot.f,aes(Cells_AreaShape_Solidity,
#                        Cells_AreaShape_Compactness,colour=as.factor(Class)))+geom_point()
# # boxplot(cell_plot.f$ALPCyPositivelasm_Intensity_IntegratedIntensity_ALP4Corr~cell_plot.f$Class)
# # boxplot(cell_plot.f$ALP_Int_scale~cell_plot.f$Class)
# 
# ggplot(cell_plot.f,aes(Cells_AreaShape_Solidity,
#                        Cells_AreaShape_Compactness,colour=ALPCyPositivelasm_Intensity_IntegratedIntensity_ALP4Corr))+geom_point()+
#   scale_colour_gradient2(limits=c(0, 1500),low="yellow", high="blue")
# 
# min(cell_plot.f$ALP_Int_scale)
# ggplot(cell_plot.f,aes(Cells_AreaShape_Solidity,
#                        Cells_AreaShape_Compactness,colour=ALP_Int_scale))+geom_point()+
#   scale_colour_gradient2(limits=c(-0.9, 4),low="yellow", high="blue")
# 
# ggplot(cell_plot.f,aes(Cells_AreaShape_Solidity,
#                        Cells_AreaShape_Compactness,colour=Class))+geom_point()
# 
#making per repaet plot()
library(plyr)
cell_plot.rf<-ddply(cell_plot.f,"ImageNumber",numcolwise(median))
cell_plot.rf.cl<-merge(cell_plot.rf,cell_plot.f[,c("ImageNumber","Class")])

# summary(cell_plot.f[cell_plot.rf.cl$Class=="pos","ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr"])
# 
# summary(cell_plot.f[cell_plot.rf.cl$Class=="neg","ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr"])

ggplot(cell_plot.rf.cl,aes(Compactness,Solidity,
         colour=as.factor(Class)))+geom_point()+theme_bw()+
  scale_color_discrete(name ="Surface Class for \n Alp expression", 
                         labels=c("Negative","Positive"))+ theme(legend.position = c(.8, .8))
# ggplot(cell_plot.rf,aes(Solidity,Compactness,colour=ALPCytoplasm_Intensity_IntegratedIntensity_ALP4Corr))+
#   geom_point()+
#   scale_color_gradient2(limits=c(0, 1000),low = "blue", midpoint = 300, mid = "yellow", high = "red",space="Lab")+
#   theme_bw()+theme(legend.position = c(.9, .9))
#  
# scale_colour_gradientn(colours=c("red","yellow","green"),
#       values  = c(0, 50, 1000))
#   
# scale_colour_gradient2(limits=c(0, 1000),low="yellow", high="blue",space="Lab")
# 
# ggplot(cell_plot.rf,aes(Compactness,
#                        Solidity,colour=ALP_Int_scale))+geom_point()+
#   scale_colour_gradient2(limits=c(0, 2),low="red",mid = "yellow", high="blue")+theme_bw()+
#   scale_color_continuous(name ="Surface Class for \n Alp expression", 
#                        labels=c("Negative","Positive"))+ theme(legend.position = c(.8, .8))
# 
# 
# 
