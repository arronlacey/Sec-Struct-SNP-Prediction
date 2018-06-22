
library(e1071)

#1st fold

svm1<-svm(label ~ ., data = f.train)

svm.table1<-table(f.test$label, predict(svm1, f.test[,1:ncol(f.test)-1], type = "class"))


#2nd fold

svm2<-svm(label ~ ., data = f2.train)

svm.table2<-table(f2.test$label, predict(svm2, f2.test[,2:ncol(f2.test)], type = "class"))


#3rd fold

svm3<-svm(label ~ ., data = f3.train)

svm.table3<-table(f3.test$label, predict(svm3, f3.test[,2:ncol(f3.test)], type = "class"))


#4th fold

svm4<-svm(label ~ ., data = f4.train)

svm.table4<-table(f4.test$label, predict(svm4, f4.test[,2:ncol(f4.test)], type = "class"))


#5th fold

svm5<-svm(label ~ ., data = f5.train)

svm.table5<-table(f5.test$label, predict(svm5, f5.test[,2:ncol(f5.test)], type = "class"))

svm.prob<-predict(svm1, f.test[,2:ncol(f2.test)],decision.values=TRUE)
svm.prob1<-attr(svm.prob,"decision.values")
svm.perf<-ROCR::performance(ROCR::prediction(1-svm.prob1,f2.test$label),'tpr','fpr')
