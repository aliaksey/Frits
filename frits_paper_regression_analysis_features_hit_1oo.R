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

###get rid of highly correlated features
corMatips <- cor(alp_model[,-c(1:6)], method="spearman")
corrplot(corMatips)

hCorr <- findCorrelation(corMatips, 0.75)
alp_model_f<-alp_model[,!colnames(alp_model)%in%colnames(corMatips[,hCorr])]
corMatips.f <- cor(alp_model_f[,-c(1:6)],method="spearman")
corrplot(corMatips.f)
##rename rows
alp_model_f[alp_model_f$feature.idx==2177,"feature.idx"]<-paste(2177,c(1:2),sep="_")
rownames(alp_model_f)<-alp_model_f$feature.idx
#Selecting hit 100 from both sides
alp_model_f.s<-alp_model_f[order(alp_model_f$median),]
alp_model_f.ss<-alp_model_f.s[c(1:100,2079:2178),]
alp_model_f.ss$Class<-c(rep("bottom",100),rep("top",100))
alp_model_f.sss<-alp_model_f.ss[,-c(1:6)]
alp_model_f.sss$Class<-as.factor(alp_model_f.sss$Class)
##creating model

##selecting samples for training and testing
data_for_model<-alp_model_f.sss
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

rpart_training2 <- rpart(Class~.,  method="class", data=alp_model_f.sss)
#plot(rpart_training)
#text(rpart_training)
##polt as party object
rpart1a2 <- as.party(rpart_training2)
plot(rpart1a2, main="Pruned CART classyfication tree for OCT4 hits")


rpartPred <- predict(rpart_training, forTesting, type = "class")
confusionMatrix(rpartPred, forTesting$Class)
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

rpartProbs <- predict(rpartTune, forTesting, type = "prob")
head(rpartProbs)
library(pROC)
rpartROC <- roc(forTesting$Class, rpartProbs[, "top"], 
                levels = rev(forTesting$Class))
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
plot(varImp(svmTune),top=17,cex=4,pch=16,cex.axis=30,
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
##compare 
plot(rpartROC, type = "S",col = "red", main="comparison 3 fitted models")#,print.thres = .5
plot(svmROC, add = TRUE, col = "blue")#,print.thres = .5
plot(logitROC, add = TRUE, col = "yellow")
legend(x=0.6,y=0.2, legend=c("CART: AUC=0.5","SVM: AUC=0.9","Logit: AUC=0.78"), 
       lty=1, col=c("red","blue","yellow"), bty='o', cex=1)
