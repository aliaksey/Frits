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
                                         !grepl("Zernike",colnames(alp_model_morph.temp))&
                                           !grepl("Orientation",colnames(alp_model_morph.temp))&
                                        !grepl("_Center_",colnames(alp_model_morph.temp))])]
colnames(alp_model_morph)

# ###get rid of highly correlated features
# corMatips <- cor(alp_model_morph[,-c(1:4)], method="spearman")
# corrplot(corMatips)
# 
# hCorr <- findCorrelation(corMatips, 0.75)
# alp_model_f<-alp_model_morph[,!colnames(alp_model_morph)%in%colnames(corMatips[,hCorr])]
# 
# alp_model_f<-alp_model_morph
# 
# colnames(alp_model_f)
# 
# corMatips.f <- cor(alp_model_f[,-c(1:4)],method="spearman")
# corrplot(corMatips.f)
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
# top.feat<-tail(unique(alp_model_f.s$feature.idx),100)
# bottom.feat<-head(unique(alp_model_f.s$feature.idx),100)
# 
# alp_model_f.ss.t<-alp_model_f.s[alp_model_f.s$feature.idx%in%top.feat,]
# alp_model_f.ss.t$Class<-"Top"
# alp_model_f.ss.b<-alp_model_f.s[alp_model_f.s$feature.idx%in%bottom.feat,]
# alp_model_f.ss.b$Class<-"Bottom"
# alp_model_f.ss<-rbind(alp_model_f.ss.b,alp_model_f.ss.t)
# plot(alp_model_f.ss.t$median)
# plot(alp_model_f.ss.b$median)
# plot(alp_model_f.ss$median)
##selecting all surfaces for analysis
alp_model_f.ss<-alp_model_f.s
##find most severe outliers (boxplot rule)
#alp_model_f.sso<-c()
#library(plyr)
rm(list=c("rslt","rsltbb","rslta"))
library(doParallel)
library(foreach)
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

rslt<-foreach(i=unique(alp_model_f.ss[,"feature.idx"]), .combine='rbind') %dopar% {
  temp<-alp_model_f.ss[alp_model_f.ss$feature.idx==i,]
  ###filter based on all features cell number
  for(f in c(5:(length(temp)))){
    temp2<-temp[,c(3,4,f)]
    #temp2$CellIdx<-paste(temp2$ImageNumber,temp2$ObjectNumber,sep="_")
    lbnda<-as.numeric(quantile(temp2[,3], probs = 0.25))
    ubnda<-as.numeric(quantile(temp2[,3], probs = 0.75))
    iuda<-ubnda-lbnda
    rslta<-temp2[temp2[,3]<(ubnda+1.5*iuda)&
                   temp2[,3]>(lbnda-1.5*iuda),]
    if (length(rslta[,1])==0) break
    if(f==5) rsltbb<-rslta else rsltbb<-merge(rslta, rsltbb, 
                                                          by=c("ImageNumber","ObjectNumber"))
   }
  rsltbb
  }
stopCluster(cl)
registerDoSEQ()
save(rslt,file="result of outliers elemination.RDATA")

alp_model_f.sso<-merge(alp_model_f.ss[,c(1:3,length(alp_model_f.ss))],rslt,by="ImageNumber") 
#alp_model_f.sso$Class<-as.factor(alp_model_f.sso$Class)
##take median per repeat
library(plyr)
alp_model_f.sso2<-ddply(alp_model_f.sso,"ImageNumber", numcolwise(median))
# alp_model_f.sso3<-merge(alp_model_f.sso2,alp_model_f.sso[,c(1,4)],by="ImageNumber")
#median per feature
alp_model_f.sso2.f<-ddply(alp_model_f.sso2,"feature.idx", numcolwise(median))

# alp_model_f.sso3.f<-merge(alp_model_f.sso2.f,unique(alp_model_f.sso[,c(2,4)]),
#                           by="feature.idx",all=F)

alp_model_f.sss.t<-alp_model_f.sso2.f[,-c(2,3,5)]
rownames(alp_model_f.sss.t)<-alp_model_f.sss.t$feature.idx
alp_model_f.sss<-alp_model_f.sss.t[,-1]
save(alp_model_f.sss,file="Median cell shape features per surface.RDATA")
################create classyfication of cells
#SCALING

cntr<-apply(alp_model_f.sss,2,function(x) median(x))
scl<-apply(alp_model_f.sss,2,function(x) mad(x))
alp_model_f.sss.scale<- scale(alp_model_f.sss,
                                 center=cntr,scale=scl)
boxplot(alp_model_f.sss.scale)
##Calculate PCA
alp_model_f.sss.pca<-prcomp(alp_model_f.sss.scale, center=F, sclale=F )

biplot(alp_model_f.sss.pca)
summary(alp_model_f.sss.pca)
plot(alp_model_f.sss.pca, type = "l")
plot(alp_model_f.sss.pca$sdev)

plot(cumsum(alp_model_f.sss.pca$sdev^2 / sum(alp_model_f.sss.pca$sdev^2)))
barplot(cumsum(alp_model_f.sss.pca$sdev^2 / sum(alp_model_f.sss.pca$sdev^2)))
barplot(alp_model_f.sss.pca$sdev^2 / sum(alp_model_f.sss.pca$sdev^2))
plot(alp_model_f.sss.pca$sdev^2 / sum(alp_model_f.sss.pca$sdev^2),type="l")
plot(alp_model_f.sss.pca$sdev^2 / sum(alp_model_f.sss.pca$sdev^2))

plot(alp_model_f.sss.pca$x[,1],alp_model_f.sss.pca$x[,2])
#calsulate distance matrix
#to.dist.cl<-alp_model_f.sss.pca$x[,1:6]
to.dist.cl<-alp_model_f.sss.scale
data.dist<-dist(to.dist.cl, method="euclidean")
#performing clustering
hclustres<-hclust(data.dist, method = "ward.D2")
plot(hclustres)
clust.numb<-4 ##specify number of clusters for selection
rect.hclust(hclustres,k=clust.numb)
library(dendextend)
hclustres.dend <- as.dendrogram(hclustres)

hclustres.dend <- color_branches(hclustres.dend, k = clust.numb)
hclustres.dend <- color_labels(hclustres.dend, k = clust.numb)
plot(hclustres.dend)

###############start catting
clstrs<-as.data.frame(cbind(Cluster=cutree(hclustres.dend,k=clust.numb,
                                           order_clusters_as_data =F,
                                           sort_cluster_numbers = T),
                            FeatureIdx=rownames(as.matrix(cutree(hclustres.dend,k=clust.numb,
                                         order_clusters_as_data =F,
                                         sort_cluster_numbers = F)))))

#clstrs<-cutree(hclustres,k=clust.numb)
##merge clusters with cell shape data
alp.model.cell.clust.t<-merge(alp_model[,-c(2:6)],clstrs,by.y="FeatureIdx",by.x="feature.idx",sort=F)
rownames(alp.model.cell.clust.t)<-alp.model.cell.clust.t$feature.idx
alp.model.cell.clust<-alp.model.cell.clust.t[,-1]
alp.model.cell.clust$Cluster<-as.factor(alp.model.cell.clust$Cluster)
levels(alp.model.cell.clust$Cluster)<-paste("Cluster",levels(alp.model.cell.clust$Cluster), sep="_")
summary(alp.model.cell.clust$Cluster)
##selecting samples for training and testing
data_for_model<-alp.model.cell.clust
##partition basd on cluster
class_data<-data_for_model[,"Cluster"]
inTrain <- createDataPartition(class_data, p =3/4, list = FALSE)
forTraining <- data_for_model[inTrain,]
#forTraining <- ipshitsmodel
forTrainingX <- forTraining[, names(forTraining) != "Cluster"]
#reate testing set for features selection
forTesting <- data_for_model[-inTrain,]

#############rpart analysis##########
#     rpart_training <- rpart(Class~.,  method="class", data=data_for_model)
#     #plot(rpart_training)
#     #text(rpart_training)
#     ##polt as party object
#     rpart1a <- as.party(rpart_training)
#     plot(rpart1a, main="Pruned CART classyfication tree for OCT4 hits")
#     
##rpart on unchanged data

rpart_training2 <- rpart(Cluster~.,  method="class", data=data_for_model)
#plot(rpart_training)
#text(rpart_training)
##polt as party object
rpart1a2 <- as.party(rpart_training2)
plot(rpart1a2, main="Pruned CART Clusteryfication tree for OCT4 hits")


# rpartPred <- predict(rpart_training, forTesting, type = "Cluster")
# confusionMatrix(rpartPred, forTesting$Cluster)
______________________________________________________________________________________ 
#Multi-Class Summary Function
#Based on caret:::twoClassSummary
require(compiler)
multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL){
  
  #Load Libraries
  require(Metrics)
  require(caret)
  
  #Check data
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  
  #Calculate custom one-vs-all stats for each class
  prob_stats <- lapply(levels(data[, "pred"]), function(class){
    
    #Grab one-vs-all data for the class
    pred <- ifelse(data[, "pred"] == class, 1, 0)
    obs  <- ifelse(data[,  "obs"] == class, 1, 0)
    prob <- data[,class]
    
    #Calculate one-vs-all AUC and logLoss and return
    cap_prob <- pmin(pmax(prob, .000001), .999999)
    prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
    names(prob_stats) <- c('ROC', 'logLoss')
    return(prob_stats) 
  })
  prob_stats <- do.call(rbind, prob_stats)
  rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"])
  
  #Aggregate and average class-wise stats
  #Todo: add weights
  class_stats <- cbind(CM$byClass, prob_stats)
  class_stats <- colMeans(class_stats)
  
  #Aggregate overall stats
  overall_stats <- c(CM$overall)
  
  #Combine overall with class-wise stats and remove some stats we don't want 
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', 
                                       'Prevalence', 'Detection Prevalence')]
  
  #Clean names and return
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  return(stats)
  
})

##tunning the model

cvCtrl <- trainControl(method = "repeatedcv", repeats = 10,
                       classProbs = TRUE,savePred=T,returnResamp="final")#,
                      #  summaryFunction = multiClassSummary)
 cl <- makeCluster(detectCores(), type='PSOCK')
 registerDoParallel(cl)

#rpart
rpartTune <- train(Cluster ~ ., data = forTraining, method = "rpart",
                   tuneLength = 10,
                   metric = 'Accuracy',
                   trControl = cvCtrl)
 stopCluster(cl)
 registerDoSEQ()
# for(stat in c('Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyPValue', 
#               'Sensitivity', 'Specificity', 'Pos_Pred_Value', 
#               'Neg_Pred_Value', 'Detection_Rate', 'ROC', 'logLoss')) {
#   
#   print(plot(rpartTune, metric=stat))
# }

plot(rpartTune)
predictors(rpartTune)
plot(varImp(rpartTune),top=5,cex=4,pch=16,
     main="Feature importance for CART method")
varImp_re<-varImp(rpartTune)
row.names(varImp_re$importance)[varImp_re$importance>0]
plot.train(rpartTune)
plot(rpartTune, scales = list(x = list(log = 10)))

rpartPred2 <- predict(rpartTune, forTesting)
confusionMatrix(rpartPred2, forTesting$Cluster)
table(forTraining$Cluster)
table(forTesting$Cluster)
##knn
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

knnTune <- train(Cluster ~ ., data = forTraining, method = "knn",
                   tuneLength = 10,
                   metric = 'Accuracy',
                   trControl = cvCtrl)
stopCluster(cl)
registerDoSEQ()
# for(stat in c('Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyPValue', 
#               'Sensitivity', 'Specificity', 'Pos_Pred_Value', 
#               'Neg_Pred_Value', 'Detection_Rate', 'ROC', 'logLoss')) {
#   
#   print(plot(knnTune, metric=stat))
# }

plot(knnTune)
predictors(knnTune)
plot(varImp(knnTune),top=10,cex=4,pch=16,
     main="Feature importance for CART method")
varImp_re<-varImp(knnTune)
plot.train(knnTune)
plot(knnTune, scales = list(x = list(log = 10)))

knnPred2 <- predict(knnTune, forTesting)
confusionMatrix(knnPred2, forTesting$Cluster)

##svm
#Setup parallel cluster
library(doParallel)
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
svmTune <- train(Cluster ~ ., data = forTraining, method = "svmRadial",
                 tuneLength = 10,
                 preProc = c("center", "scale"),
                   metric = 'Accuracy',
                   trControl = cvCtrl)
# for(stat in c('Accuracy', 'Kappa', 'AccuracyLower', 'AccuracyUpper', 'AccuracyPValue', 
#               'Sensitivity', 'Specificity', 'Pos_Pred_Value', 
#               'Neg_Pred_Value', 'Detection_Rate', 'ROC', 'logLoss')) {
#   
#   print(plot(svmTune, metric=stat))
# }
stopCluster(cl)
registerDoSEQ()

plot(svmTune)
predictors(svmTune)
plot(varImp(svmTune),top=10,cex=4,pch=16,
     main="Feature importance for CART method")
varImp_re<-varImp(svmTune)
row.names(varImp_re$importance)[varImp_re$importance>0]
plot.train(svmTune)
plot(svmTune, scales = list(x = list(log = 10)))

svmPred2 <- predict(svmTune, forTesting)
confusionMatrix(svmPred2, forTesting$Cluster)

##logit
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

logitTune <- train(x = forTrainingX,
                   y = forTraining$Cluster,
                   method = "glm",
                   tuneLength = 10,
                   family = binomial(link = "logit"),
                   metric = 'Accuracy',
                   trControl = cvCtrl)
stopCluster(cl)
registerDoSEQ()
summary(logitTune)

#predictors(logitTune)
plot(varImp(logitTune),top=10,cex=4,pch=16,cex.axis=30,
     main="Feature importance for Logit method")
logitTune$finalModel
logitPred <- predict(logitTune, forTesting[, names(forTesting) != "Cluster"])
confusionMatrix(logitPred, forTesting$Cluster)

##nb
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
nbTune <- train(x = forTrainingX,
                y = forTraining$Cluster,
                method = "nb",
                tuneLength = 10,
                metric = 'Accuracy',
                trControl = cvCtrl)

stopCluster(cl)
registerDoSEQ()
summary(nbTune)

predictors(nbTune)
plot(varImp(nbTune),
     main="Feature importance for NB method")
nbTune$finalModel
nbPred <- predict(nbTune, forTesting[, names(forTesting) != "Cluster"])
confusionMatrix(nbPred, forTesting$Cluster)

##rf
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
rfTune <- train(x = forTrainingX,
                y = forTraining$Cluster,
                method = "rf",
                tuneLength = 10,
                allowParallel=TRUE,
                metric = 'Accuracy',
                trControl = cvCtrl)

stopCluster(cl)
registerDoSEQ()

summary(rfTune)

predictors(rfTune)
plot(varImp(rfTune),
     main="Feature importance for Logit method")
rfTune$finalModel
rfPred <- predict(rfTune, forTesting[, names(forTesting) != "Cluster"])
confusionMatrix(rfPred, forTesting$Cluster)

##LDA
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

ldaTune <- train(x = forTrainingX,
                 y = forTraining$Cluster,
                 method = "lda",
                 tuneLength = 10,
                 metric = 'Accuracy',
                 trControl = cvCtrl)
stopCluster(cl)
registerDoSEQ()
summary(ldaTune)

predictors(ldaTune)
plot(varImp(ldaTune),top=10,cex=4,pch=16,cex.axis=30,
     main="Feature importance for lda method")
ldaTune$finalModel
ldaPred <- predict(ldaTune, forTesting[, names(forTesting) != "Cluster"])
confusionMatrix(ldaPred, forTesting$Cluster)

##QDA
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

qdaTune <- train(x = forTrainingX,
                 y = forTraining$Cluster,
                 method = "qda",
                 tuneLength = 10,
                 metric = 'Accuracy',
                 trControl = cvCtrl)


stopCluster(cl)
registerDoSEQ()
summary(qdaTune)

predictors(qdaTune)
plot(varImp(qdaTune),top=10,cex=4,pch=16,cex.axis=30,
     main="Feature importance for qda method")
qdaTune$finalModel
qdaPred <- predict(qdaTune, forTesting[, names(forTesting) != "Cluster"])
confusionMatrix(qdaPred, forTesting$Cluster)
#########################################################################################
##########################################################################################
##########################################################################################
############################################################################################
#################checking regression model does it works better
#alp_model_f.sss.t
##selecting feature to use as predictor
selected_variable<-"Cells_AreaShape_FormFactor"
alp.model.cell.regr.t<-merge(alp_model[,-c(2:6)],alp_model_f.sss.t[,
           c("feature.idx",selected_variable)],by="feature.idx",sort=F)
colnames(alp.model.cell.regr.t)[colnames(alp.model.cell.regr.t)==selected_variable]<-
  "Cell_Shape_Feature"
alp.model.cell.regr<-alp.model.cell.regr.t[,-1]

##selecting samples for training and testing
data_for_model<-alp.model.cell.regr

class_data<-data_for_model$Cell_Shape_Feature
inTrain <- createDataPartition(class_data, p =3/4, list = FALSE)
forTraining <- data_for_model[inTrain,]
forTrainingX <- forTraining[, !colnames(forTraining)=="Cell_Shape_Feature"]
#create testing set for features selection
forTesting <- data_for_model[-inTrain,]

##regression tree with CART metod for whole data

rpart_training2 <-rpart(Cell_Shape_Feature~.,  method="anova", 
      data=data_for_model)
# plot(rpart_training2)
# text(rpart_training2)
##polt as party object
rpart1a2 <- as.party(rpart_training2)
plot(rpart1a2, main="anova with CART for Alp integrated intensity")

##regression tree with CART metod for simplified data
rpart_training3 <- rpart(Cell_Shape_Feature~.,  method="anova", 
                         data=data_for_model)
printcp(rpart_training3)
plotcp(rpart_training3)
rsq.rpart(rpart_training3)
summary(rpart_training3)
rpart_training4<-prune(rpart_training3, cp= rpart_training3$cptable[which.min(rpart_training3$cptable[,"xerror"]),"CP"]) 
##polt as party object
rpart1a3 <- as.party(rpart_training3)
plot(rpart1a3, main="anova with CART for Alp integrated intensity")
rpart1a4 <- as.party(rpart_training4)
plot(rpart1a4, main="anova with CART for Alp integrated intensity,Pruned")


# rpartPred <- predict(rpart_training3, forTesting)
# summary(rpartPred)
#______________________________________________________________________________________ 
##tunning the model
cvCtrl <- trainControl(method = "repeatedcv", repeats = 10,
                       summaryFunction = defaultSummary,
                       savePred=T,returnResamp="final")
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
#rpart
rpartTune <- train(Cell_Shape_Feature ~ ., data = forTraining, method = "rpart",
                   tuneLength = 10,
                   metric = "RMSE",
                   trControl = cvCtrl)
stopCluster(cl)
registerDoSEQ()
plot(rpartTune)
predictors(rpartTune)
plot(varImp(rpartTune),top=10,cex=4,pch=16,
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
Observed = forTesting$Cell_Shape_Feature
Predicted = predict(rpartTune, forTesting)
library(latticeExtra)
xyplot(Observed ~ Predicted, panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.lmlineq(x, y, adj = c(1,0), lty = 1,xol.text='red',
                col.line = "blue", digits = 1,r.squared =TRUE)
})

library(doParallel)
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
##svm
svmTune <- train(Cell_Shape_Feature ~ ., data = forTraining,
                 method = "svmRadial",
                 tuneLength = 10,
                 preProc = c("center", "scale"),
                 metric = "RMSE",
                 trControl = cvCtrl)
stopCluster(cl)
registerDoSEQ()

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
Observed = forTesting$Cell_Shape_Feature
Predicted = predict(svmTune, forTesting)
library(latticeExtra)
xyplot(Observed ~ Predicted, panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.lmlineq(x, y, adj = c(1,0), lty = 1,xol.text='red',
                col.line = "blue", digits = 1,r.squared =TRUE)
})
##RF
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
rfTune <- train(x = forTrainingX,
                y = forTraining$Cell_Shape_Feature,
                method = "rf",
                tuneLength = 10,
                prox=TRUE,
                allowParallel=TRUE,
                metric = "RMSE",
                trControl = cvCtrl)
stopCluster(cl)
registerDoSEQ()
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
Observed = forTesting$Cell_Shape_Feature
Predicted = predict(rfTune, forTesting)
library(latticeExtra)
cor(Observed,Predicted)
xyplot(Observed ~ Predicted, 
       panel = function(x, y, ...) {
  panel.xyplot(x, y, ...)
  panel.lmlineq(x, y, adj = c(1,0), lty = 1,xol.text='red',
                col.line = "blue", digits = 1,r.squared =TRUE)
})

##compare 
predTargets.c <- extractPrediction(list(svmTune,rpartTune,rfTune), forTesting)
plotObsVsPred(predTargets.c)
