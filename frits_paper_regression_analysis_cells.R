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
##loading morphological data from phenome 
load("D:/projects/phenome_project/phenome/Cell all data & ground truth scaled.RData")
#attach alp rank to morphological data
alp_model_morph.temp<-merge(alp_model,feature.cell.scale,by.x="feature.idx",by.y="FeatureIdx")
#Selecting features for anlaysis
alp_model_morph<-alp_model_morph.temp[,c("feature.idx","median",
              colnames(alp_model_morph.temp)[grepl("Cells_AreaShape_",colnames(alp_model_morph.temp))&
                                               !grepl("_Center_",colnames(alp_model_morph.temp))])]
colnames(alp_model_morph)
# ###get rid of highly correlated features
# corMatips <- cor(alp_model[,-c(1:6)], method="spearman")
# corrplot(corMatips)
# 
# hCorr <- findCorrelation(corMatips, 0.75)
# alp_model_f<-alp_model[,!colnames(alp_model)%in%colnames(corMatips[,hCorr])]
# corMatips.f <- cor(alp_model_f[,-c(1:6)],method="spearman")
# corrplot(corMatips.f)
##rename rows
alp_model_morph[alp_model_morph$feature.idx==2177,"feature.idx"]<-paste(2177,c(1:2),sep="_")
rownames(alp_model_morph)<-alp_model_morph$feature.idx
##creating model

##selecting samples for training and testing
data_for_model<-alp_model_morph
variables_class<-c("feature.idx", "ranksum","median","mean","ttest", "mestimation" )

class_data<-data_for_model[,variables_class[3]]
inTrain <- createDataPartition(class_data, p =3/4, list = FALSE)
forTraining <- data_for_model[inTrain,!colnames(data_for_model)%in%variables_class[-3]]
forTrainingX <- forTraining[, !colnames(forTraining)%in%variables_class[3]]
#create testing set for features selection
forTesting <- data_for_model[-inTrain,!colnames(data_for_model)%in%variables_class[-3]]

##regression tree with CART metod for whole data
as.symbol(variables_class[3])
rpart_training2 <-                                    rpart(median~.,  method="anova", 
                                                            data=as.data.frame(data_for_model[,!colnames(data_for_model)%in%variables_class[-3]]))
# plot(rpart_training2)
# text(rpart_training2)
##polt as party object
rpart1a2 <- as.party(rpart_training2)
plot(rpart1a2, main="anova with CART for Alp integrated intensity")
______________________________________________________________________________________ 
##tunning the model
cvCtrl <- trainControl(method = "repeatedcv", repeats = 10,
                       summaryFunction = defaultSummary,
                       savePred=T,returnResamp="final")
#rpart
rpartTune <- train(median ~ ., data = forTraining, method = "rpart",
                   tuneLength = 10,
                   metric = "RMSE",
                   trControl = cvCtrl)
plot(rpartTune)
predictors(rpartTune)
plot(varImp(rpartTune),top=5,cex=4,pch=16,
     main="Feature importance for CART method")
varImp_re<-varImp(rpartTune)
row.names(varImp_re$importance)[varImp_re$importance>0]
plot.train(rpartTune)
#print.train(rpartTune)
plot(rpartTune, scales = list(x = list(log = 10)))
##plot prediction vs observed
#in caret package 
predTargets <- extractPrediction(list(rpartTune), testX=forTesting)
plotObsVsPred(predTargets)
#in latice extra package
Observed = forTesting$median
Predicted = predict(rpartTune, forTesting)
library(latticeExtra)
xyplot(Observed ~ Predicted, panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.lmlineq(x, y, adj = c(1,0), lty = 1,xol.text='red',
                col.line = "blue", digits = 1,r.squared =TRUE)
})

##svm
svmTune <- train(median ~ ., data = forTraining,
                 method = "svmRadial",
                 tuneLength = 10,
                 preProc = c("center", "scale"),
                 metric = "RMSE",
                 trControl = cvCtrl)
svmTune

#predictors(svmTune)
plot(varImp(svmTune),top=17,cex=4,pch=16,cex.axis=30,
     main="Feature importance for SVM method")
svmTune$finalModel
##plot prediction vs observed
#in caret package 
predTargets <- extractPrediction(list(svmTune), testX=forTesting)
plotObsVsPred(predTargets)
#in latice extra package
Observed = forTesting$median
Predicted = predict(svmTune, forTesting)
library(latticeExtra)
xyplot(Observed ~ Predicted, panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.lmlineq(x, y, adj = c(1,0), lty = 1,xol.text='red',
                col.line = "blue", digits = 1,r.squared =TRUE)
})
##RF
rfTune <- train(x = forTrainingX,
                y = forTraining$median,
                method = "rf",
                tuneLength = 10,
                prox=TRUE,
                allowParallel=TRUE,
                metric = "RMSE",
                trControl = cvCtrl)
summary(rfTune)


#predictors(rfTune)
plot(varImp(rfTune),top=17,cex=4,pch=16,cex.axis=30,
     main="Feature importance for rf method")
rfTune$finalModel
##plot prediction vs observed
#in caret package 
predTargets <- extractPrediction(list(rfTune), testX=forTesting)
plotObsVsPred(predTargets)
#in latice extra package
Observed = forTesting$median
Predicted = predict(rfTune, forTesting)
library(latticeExtra)
xyplot(Observed ~ Predicted, panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.lmlineq(x, y, adj = c(1,0), lty = 1,xol.text='red',
                col.line = "blue", digits = 1,r.squared =TRUE)
})

##compare 
predTargets.c <- extractPrediction(list(svmTune,rpartTune,rfTune), forTesting)
plotObsVsPred(predTargets.c)
