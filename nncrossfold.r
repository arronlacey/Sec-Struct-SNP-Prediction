
library(nnet)

#1st fold

nn<-nnet(label ~ ., data = f.train, size = 5, rang = 0.1,
         decay = 5e-4, maxit = 500)

nn.table1<-table(f.test$label, predict(nn, f.test[,2:ncol(f.test)], type = "class"))


#2nd fold

nn2<-nnet(label ~ ., data = f2.train, size = 5, rang = 0.1,
         decay = 5e-4, maxit = 500)

nn.table2<-table(f2.test$label, predict(nn2, f2.test[,2:ncol(f2.test)], type = "class"))


#3rd fold

nn3<-nnet(label ~ ., data = f3.train, size = 5, rang = 0.1,
         decay = 5e-4, maxit = 500)

nn.table3<-table(f3.test$label, predict(nn3, f3.test[,2:ncol(f3.test)], type = "class"))


#4th fold

nn4<-nnet(label ~ ., data = f4.train, size = 5, rang = 0.1,
         decay = 5e-4, maxit = 500)

nn.table4<-table(f4.test$label, predict(nn4, f4.test[,2:ncol(f4.test)], type = "class"))


#5th fold

nn5<-nnet(label ~ ., data = f5.train, size = 5, rang = 0.1,
         decay = 5e-4, maxit = 500)

nn.table5<-table(f5.test$label, predict(nn5, f5.test[,2:ncol(f5.test)], type = "class"))


nn.prob<-predict(nn, f.test,type = "raw")
nn.perf<-ROCR::performance(prediction(nn.prob,f.test$label),'tpr','fpr')