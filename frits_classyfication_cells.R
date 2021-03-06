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
alp_model_f<-alp_model_morph
# 
# 
# 
# ###get rid of highly correlated features
# corMatips <- cor(alp_model[,-c(1:6)], method="spearman")
# corrplot(corMatips)
# 
# hCorr <- findCorrelation(corMatips, 0.75)
# alp_model_f<-alp_model[,!colnames(alp_model)%in%colnames(corMatips[,hCorr])]
# corMatips.f <- cor(alp_model_f[,-c(1:6)],method="spearman")
# corrplot(corMatips.f)
##rename rows
alp_model_f[alp_model_f$feature.idx==2177,"feature.idx"]<-paste(2177,c(1:2),sep="_")
rownames(alp_model_f)<-alp_model_f$feature.idx
#add new ranking
load("Alp_intensities_rank_Alex.RData")
alp_model_f<-merge(alp_model_f,alp_integr_intensities.f,by.x="feature.idx",by.y="FeatureIdx")
alp_model_f$median.new<-alp_model_f$AlpTrMeanFeatureMean
alp_model_f<-alp_model_f[,!colnames(alp_model_f)%in%c("AlpTrMeanFeatureMean","AlpTrMeanFeatureMedian")]
#Selecting hit 100 from both sides
alp_model_f.s<-alp_model_f[order(alp_model_f$median),]
alp_model_f.ss<-alp_model_f.s[c(1:100,2071:2170),]
alp_model_f.ss$Class<-c(rep("bottom",100),rep("top",100))
alp_model_f.sss<-alp_model_f.ss[,-c(1:6)]
alp_model_f.sss$Class<-as.factor(alp_model_f.sss$Class)
#for new calculation
alp_model_f.s2<-alp_model_f[order(alp_model_f$median.new),]
alp_model_f.ss2<-alp_model_f.s2[c(1:100,2071:2170),]
alp_model_f.ss2$Class<-c(rep("bottom",100),rep("top",100))
alp_model_f.sss2<-alp_model_f.ss2[,-c(1:6)]
alp_model_f.sss2$Class<-as.factor(alp_model_f.sss2$Class)
##check are any correlation between two ranks
plot(alp_model_f.sss$median,alp_model_f.sss2$median.new)
table(alp_model_f.sss2$median.new)
# ##selecting hits
# x<-row.names(alp_model_f.sss[alp_model_f.sss$Class=="top",])
# writeClipboard(x)
# save(x, file="alp.high.RDATA")
# 
# y<-row.names(alp_model_f.sss[alp_model_f.sss$Class=="bottom",])
# 
# save(y, file="alp.low.RDATA")
# 

##creating model

##selecting samples for training and testing
data_for_model<-alp_model_f.sss[,colnames(alp_model_f.sss)!="median.new"]
#data_for_model<-alp_model_f.sss[,colnames(alp_model_f.sss2)!="median"]
class_data<-data_for_model[,"Class"]
inTrain <- createDataPartition(class_data, p =3/4, list = FALSE)
forTraining <- data_for_model[inTrain,]
#forTraining <- ipshitsmodel
forTrainingX <- forTraining[, names(forTraining) != "Class"]
#reate testing set for features selection
forTesting <- data_for_model[-inTrain,]
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
plot(rpart1a2, main="Pruned CART classyfication tree for OCT4 hits")


# rpartPred <- predict(rpart_training, forTesting, type = "class")
# confusionMatrix(rpartPred, forTesting$Class)
______________________________________________________________________________________ 
##tunning the model
cvCtrl <- trainControl(method = "repeatedcv", repeats = 10,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,savePred=T,returnResamp="final")
#rpart
rpartTune <- train(Class ~ ., data = forTraining, method = "rpart",
                   tuneLength = 10,
                   metric = "ROC",
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

rpartPred2 <- predict(rpartTune, forTesting)
confusionMatrix(rpartPred2, forTesting$Class)
table(forTesting$Class)

rpartProbs <- predict(rpartTune, forTesting, type = "prob")
head(rpartProbs)
library(pROC)
rpartROC <- roc(predictor = rpartProbs$top,
                response = forTesting$Class,
                levels = rev(levels(forTesting$Class)))

plot(rpartROC, type = "S", print.thres = .5)
rpartROC
##svm
svmTune <- train(Class ~ ., data = forTraining,
                 method = "svmRadial",
                 tuneLength = 10,
                 preProc = c("center", "scale"),
                 metric = "ROC",
                 trControl = cvCtrl)#,
#classProbs =  TRUE)
svmTune

#predictors(svmTune)
plot(varImp(svmTune),
     main="Feature importance for SVM method")
svmTune$finalModel
plot(svmTune, metric = "ROC", scales = list(x = list(log =2)))
svmPred <- predict(svmTune, forTesting[, names(forTesting) != "Class"])
confusionMatrix(svmPred, forTesting$Class)

svmProbs <- predict(svmTune, forTesting[, names(forTesting) != "Class"],
                    type = "prob")
head(svmProbs)
svmROC <- roc(predictor = svmProbs$top,
              response = forTesting$Class,
              levels = rev(levels(forTesting$Class)))
plot(svmROC, type = "S", print.thres = .5)
svmROC

##logit
logitTune <- train(x = forTrainingX,
                   y = forTraining$Class,
                   method = "glm",
                   tuneLength = 10,
                   family = binomial(link = "logit"),
                   metric = "ROC",
                   trControl = cvCtrl)
summary(logitTune)

#nb vs glm vs lda vs svm vs rf

#predictors(logitTune)
plot(varImp(logitTune),top=17,cex=4,pch=16,cex.axis=30,
     main="Feature importance for Logit method")
logitTune$finalModel
logitPred <- predict(logitTune, forTesting[, names(forTesting) != "Class"])
confusionMatrix(logitPred, forTesting$Class)

logitProbs <- predict(logitTune, forTesting[, names(forTesting) != "Class"],
                      type = "prob")
head(logitProbs)
logitROC <- roc(predictor = logitProbs$top,
                response = forTesting$Class,
                levels = rev(levels(forTesting$Class)))
plot(logitROC, type = "S", print.thres = .5)

##nb
nbTune <- train(x = forTrainingX,
                y = forTraining$Class,
                method = "nb",
                tuneLength = 10,
                metric = "ROC",
                trControl = cvCtrl)
summary(nbTune)

predictors(nbTune)
plot(varImp(nbTune),
     main="Feature importance for NB method")
nbTune$finalModel
nbPred <- predict(nbTune, forTesting[, names(forTesting) != "Class"])
confusionMatrix(nbPred, forTesting$Class)

nbProbs <- predict(nbTune, forTesting[, names(forTesting) != "Class"],
                   type = "prob")
head(nbProbs)
nbROC <- roc(predictor = nbProbs$top,
             response = forTesting$Class,
             levels = rev(levels(forTesting$Class)))
plot(nbROC, type = "S", print.thres = .5)

##rf
rfTune <- train(x = forTrainingX,
                y = forTraining$Class,
                method = "rf",
                tuneLength = 10,
                allowParallel=TRUE,
                metric = "ROC",
                trControl = cvCtrl)

summary(rfTune)

predictors(rfTune)
plot(varImp(rfTune),
     main="Feature importance for Logit method")
rfTune$finalModel
rfPred <- predict(rfTune, forTesting[, names(forTesting) != "Class"])
confusionMatrix(rfPred, forTesting$Class)

rfProbs <- predict(rfTune, forTesting[, names(forTesting) != "Class"],
                   type = "prob")
head(rfProbs)
rfROC <- roc(predictor = rfProbs$top,
             response = forTesting$Class,
             levels = rev(levels(forTesting$Class)))
plot(rfROC, type = "S", print.thres = .5)

##LDA
ldaTune <- train(x = forTrainingX,
                 y = forTraining$Class,
                 method = "lda",
                 tuneLength = 10,
                 metric = "ROC",
                 trControl = cvCtrl)
summary(ldaTune)

predictors(ldaTune)
plot(varImp(ldaTune),top=17,cex=4,pch=16,cex.axis=30,
     main="Feature importance for lda method")
ldaTune$finalModel
ldaPred <- predict(ldaTune, forTesting[, names(forTesting) != "Class"])
confusionMatrix(ldaPred, forTesting$Class)

ldaProbs <- predict(ldaTune, forTesting[, names(forTesting) != "Class"],
                    type = "prob")
head(ldaProbs)
ldaROC <- roc(predictor = ldaProbs$top,
              response = forTesting$Class,
              levels = rev(levels(forTesting$Class)))
plot(ldaROC, type = "S", print.thres = .5)

##QDA
qdaTune <- train(x = forTrainingX,
                 y = forTraining$Class,
                 method = "qda",
                 tuneLength = 10,
                 metric = "ROC",
                 trControl = cvCtrl)
summary(qdaTune)

predictors(qdaTune)
plot(varImp(qdaTune),top=17,cex=4,pch=16,cex.axis=30,
     main="Feature importance for qda method")
qdaTune$finalModel
qdaPred <- predict(qdaTune, forTesting[, names(forTesting) != "Class"])
confusionMatrix(qdaPred, forTesting$Class)

qdaProbs <- predict(qdaTune, forTesting[, names(forTesting) != "Class"],
                    type = "prob")
head(qdaProbs)
qdaROC <- roc(predictor = qdaProbs$top,
              response = forTesting$Class,
              levels = rev(levels(forTesting$Class)))
plot(qdaROC, type = "S", print.thres = .5)


##compare 
plot(rpartROC, type = "S",col = "red", main="comparison 3 fitted models")#,print.thres = .5
plot(svmROC, add = TRUE, col = "blue")#,print.thres = .5
plot(logitROC, add = TRUE, col = "yellow")
plot(nbROC, add = TRUE, col = "grey")
plot(rfROC, add = TRUE, col = "green")
plot(ldaROC, add = TRUE, col = "black")
plot(qdaROC, add = TRUE, col = "purple")
legend(x=0.45,y=0.5, legend=c("CART: AUC=0.72","SVM: AUC=0.88",
                              "Logit: AUC=0.79","NB: AUC=0.87",
                              "RF: AUC=0.86","LDA: AUC=0.87","QDA: AUC=0.72"), 
       lty=1, col=c("red","blue","yellow","grey","green","black","purple"), bty='o', cex=1)
###find most important parameters based on methods an
#For SVM classification models, the default behavior 
#is to compute the area under the ROC curve. 
plot(varImp(svmTune))
plot(varImp(svmTune,useModel = T,scale = FALSE))
#area under the ROC curve for each predictor
rocimpsvm<-filterVarImp(x = forTraining[, -ncol(forTraining)], y = forTraining$Class)
rocimpsvm[order(-rocimpsvm[,2]),]

## plot alp intensity values base on that parameters
library(ggplot2)
ggplot(data_for_model,aes(Cells_AreaShape_Solidity,
                          Cells_AreaShape_Zernike_2_2,colour=as.factor(Class)))+geom_point()
