
library(e1071)

#1st fold

nb<-naiveBayes(label ~ ., data = f.train)

nb.table1<-table(f.test$label, predict(nb, f.test[,2:ncol(f.test)], type = "class"))
nb.prob1<-predict(nb, f.test[,2:ncol(f.test)],type="raw")
nb.perf1<-ROCR::performance(ROCR::prediction(nb.prob1[,2],f.test$label),'tpr','fpr')

#2nd fold

nb2<-naiveBayes(label ~ ., data = f2.train)

nb.table2<-table(f2.test$label, predict(nb2, f2.test[,2:ncol(f2.test)], type = "class"))
nb.prob2<-predict(nb, f2.test[,2:ncol(f2.test)],type="raw")
nb.perf2<-ROCR::performance(ROCR::prediction(nb.prob2[,2],f2.test$label),'tpr','fpr')

#3rd fold

nb3<-naiveBayes(label ~ ., data = f3.train)

nb.table3<-table(f3.test$label, predict(nb3, f3.test[,2:ncol(f3.test)], type = "class"))
nb.prob3<-predict(nb, f3.test[,2:ncol(f3.test)],type="raw")
nb.perf3<-ROCR::performance(ROCR::prediction(nb.prob3[,2],f3.test$label),'tpr','fpr')

#4th fold

nb4<-naiveBayes(label ~ ., data = f4.train)

nb.table4<-table(f4.test$label, predict(nb4, f4.test[,2:ncol(f4.test)], type = "class"))
nb.prob4<-predict(nb, f4.test[,2:ncol(f4.test)],type="raw")
nb.perf4<-ROCR::performance(ROCR::prediction(nb.prob4[,2],f4.test$label),'tpr','fpr')

#5th fold

nb5<-naiveBayes(label ~ ., data = f5.train)

nb.table5<-table(f5.test$label, predict(nb5, f5.test[,2:ncol(f5.test)], type = "class"))

nb.prob5<-predict(nb, f5.test[,2:ncol(f5.test)],type="raw")
nb.perf5<-ROCR::performance(ROCR::prediction(nb.prob5[,2],f5.test$label),'tpr','fpr')
