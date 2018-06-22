
library(rpart)

#1st fold

dt1<-rpart(label ~ .,method="class", data = f.train)

dt.table1<-table(f.test$label, predict(dt1, f.test[,2:ncol(f.test)], type = "class"))


#2nd fold

dt2<-rpart(label ~ ., data = f2.train)

dt.table2<-table(f2.test$label, predict(dt2, f2.test[,2:ncol(f2.test)], type = "class"))


#3rd fold

dt3<-rpart(label ~ ., data = f3.train)

dt.table3<-table(f3.test$label, predict(dt3, f3.test[,2:ncol(f3.test)], type = "class"))


#4th fold

dt4<-rpart(label ~ ., data = f4.train)

dt.table4<-table(f4.test$label, predict(dt4, f4.test[,2:ncol(f4.test)], type = "class"))


#5th fold

dt5<-rpart(label ~ ., data = f5.train)

dt.table5<-table(f5.test$label, predict(dt5, f5.test[,2:ncol(f5.test)], type = "class"))


dt.prob<-predict(dt1, f.test[,2:ncol(f.test)],type="prob")
dt.perf<-ROCR::performance(ROCR::prediction(dt.prob[,2],f.test$label),'tpr','fpr')

