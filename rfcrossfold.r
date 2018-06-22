# Cross fold validation

# 1st fold

#bound <- floor((nrow(f2)/4)*3)         ORIGINAL DATASET IN THESIS
bound <- floor((nrow(df2)/10)*9)
#df <- f2[sample(nrow(f2)), ]           #shuffles the data  ORIGINAL DATASET IN THESIS
df <- df2[sample(nrow(df2)), ] 
f.train <- df[1:bound, ]              #get training set
f.test <- df[(bound+1):nrow(df), ]    #get test set

set.seed(415)

require(randomForest)
rf.mdl <- randomForest(label ~.,data=f.train,importance=TRUE)
rf.prob<-predict(rf.mdl, f.test,type = "prob")
rf.pd<-predict(rf.mdl, f.test)
rf.table1<-table(observed = f.test$label, predicted = rf.pd)
varImpPlot(rf.mdl)


# 2nd fold

#bound <- floor((nrow(f2)/4)*3)         ORIGINAL DATASET IN THESIS
bound <- floor((nrow(df2)/4)*3)
#df <- f2[sample(nrow(f2)), ]           #shuffles the data  ORIGINAL DATASET IN THESIS
df <- df2[sample(nrow(df2)), ] 
f2.train <- df[1:bound, ]              #get training set
f2.test <- df[(bound+1):nrow(df), ]    #get test set

set.seed(415)

require(randomForest)
rf.mdl2 <- randomForest(label ~.,data=f2.train,importance=TRUE)
rf.prob2<-predict(rf.mdl2, f2.test,type = "prob")
rf.pd2<-predict(rf.mdl2, f2.test)
rf.table2<-table(observed = f2.test$label, predicted = rf.pd2)






# 3rd fold

#bound <- floor((nrow(f2)/4)*3)         ORIGINAL DATASET IN THESIS
bound <- floor((nrow(df2)/4)*3)
#df <- f2[sample(nrow(f2)), ]           #shuffles the data  ORIGINAL DATASET IN THESIS
df <- df2[sample(nrow(df2)), ] 
f3.train <- df[1:bound, ]              #get training set
f3.test <- df[(bound+1):nrow(df), ]    #get test set

set.seed(415)

require(randomForest)
rf.mdl3 <- randomForest(label ~.,data=f3.train,importance=TRUE)
rf.prob3<-predict(rf.mdl3, f3.test,type = "prob")
rf.pd3<-predict(rf.mdl3, f3.test)
rf.table3<-table(observed = f3.test$label, predicted = rf.pd3)







# 4th fold

#bound <- floor((nrow(f2)/4)*3)         ORIGINAL DATASET IN THESIS
bound <- floor((nrow(df2)/4)*3)
#df <- f2[sample(nrow(f2)), ]           #shuffles the data  ORIGINAL DATASET IN THESIS
df <- df2[sample(nrow(df2)), ] 
f4.train <- df[1:bound, ]              #get training set
f4.test <- df[(bound+1):nrow(df), ]    #get test set

set.seed(415)

require(randomForest)
rf.mdl4 <- randomForest(label ~.,data=f4.train,importance=TRUE)
rf.prob4<-predict(rf.mdl4, f4.test,type = "prob")
rf.pd4<-predict(rf.mdl4, f4.test)
rf.table4<-table(observed = f4.test$label, predicted = rf.pd4)






# 5th fold

#bound <- floor((nrow(f2)/4)*3)         ORIGINAL DATASET IN THESIS
bound <- floor((nrow(df2)/4)*3)
#df <- f2[sample(nrow(f2)), ]           #shuffles the data  ORIGINAL DATASET IN THESIS
df <- df2[sample(nrow(df2)), ] 
f5.train <- df[1:bound, ]              #get training set
f5.test <- df[(bound+1):nrow(df), ]    #get test set

set.seed(415)

require(randomForest)
rf.mdl5 <- randomForest(label ~.,data=f5.train,importance=TRUE)
rf.prob5<-predict(rf.mdl5, f5.test,type = "prob")
rf.pd5<-predict(rf.mdl5, f5.test)
rf.table5<-table(observed = f5.test$label, predicted = rf.pd5)

######################## importance

rf.imp<-as.data.frame(importance(rf.mdl4,type=1,class=c("Disease","Polymorphism")))
rf.imp <- cbind(Feature = rownames(rf.imp), rf.imp)
rf.imp<-data.frame(rf.imp$Feature,rf.imp$Disease,rf.imp$Polymorphism)
names(rf.imp)<-c("Feature","Disease","Polymorphism")
rf.imp<-rf.imp[order(-rf.imp$Disease),]
rf.imp$Feature<-gsub(".1", "", rf.imp$Feature)
rf.imp[rf.imp$Feature == "diff2",1]<-"SecStruct Change Helix"
rf.imp[rf.imp$Feature == "dif",1]<-"SecStruct Change Coil"
rf.imp[rf.imp$Feature == "diff3",1]<-"SecStruct Change Sheet"
rf.imp<-rf.imp[!duplicated(rf.imp$Feature), ]
rf.imp$Feature <- factor(rf.imp$Feature, levels = rf.imp$Feature[order(rf.imp$Disease)])
rf.imp<-rf.imp[, -grep("1$", rf.imp$Feature)]


ggplot() + geom_point(data=rf.imp, mapping = aes(y=Feature,x=Disease)) +
  theme_bw() +
  ggtitle("Feature Importance Plot") +
  xlab("Mean Decrease Accuracy")
#######################################################

library(ROCR)
rf.perf<-ROCR::performance(prediction(rf.prob4[,2],f4.test$label),'tpr','fpr')
rf.auc<-ROCR::performance(ROCR::prediction(rf.prob4[,2],f4.test$label)
                          ,'tpr','fpr',measure="auc")
cutoffs <- data.frame(cut=rf.perf@alpha.values[[1]], fpr=rf.perf@x.values[[1]], 
                      tpr=rf.perf@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
stats<-cbind("Random Forest",head(subset(cutoffs, tpr <= 0.95),n=1),unlist(rf.auc@y.values))
names(stats)<-c("Classifier","Threshold","FPR","TPR","AUC")
stats$Threshold<-round(stats$Threshold,2)
stats$FPR<-round(stats$FPR,2)
stats$AUC<-round(stats$AUC,2)
stats