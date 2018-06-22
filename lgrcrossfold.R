
#1st fold

lgr<-glm(label ~ .,family = binomial("logit"), data = f.train)

lgr$score<-predict(lgr, f.test[,2:ncol(f.test)],type="response")
pred_df <- data.frame(dep_var = f.test$label, score = lgr$score)
pred_df$pred<-pred_df$score>0.5
lgr.table1<-table(pred_df$dep_var,pred_df$pred)

#2nd fold

lgr2<-glm(label ~ ., family = binomial("logit"),data = f2.train)

lgr2$score<-predict(lgr2, f2.test[,2:ncol(f2.test)],type="response")
pred_df2 <- data.frame(dep_var = f2.test$label, score = lgr2$score)
pred_df2$pred<-pred_df2$score>0.5
lgr.table2<-table(pred_df2$dep_var,pred_df2$pred)


#3rd fold

lgr3<-glm(label ~ ., family = binomial("logit"),data = f3.train)

lgr3$score<-predict(lgr3, f3.test[,2:ncol(f3.test)],type="response")
pred_df3 <- data.frame(dep_var = f3.test$label, score = lgr3$score)
pred_df3$pred<-pred_df3$score>0.5
lgr.table3<-table(pred_df3$dep_var,pred_df3$pred)


#4th fold

lgr4<-glm(label ~ .,family = binomial("logit"), data = f4.train)

lgr4$score<-predict(lgr4, f4.test[,2:ncol(f4.test)],type="response")
pred_df4 <- data.frame(dep_var = f4.test$label, score = lgr4$score)
pred_df4$pred<-pred_df4$score>0.5
lgr.table4<-table(pred_df4$dep_var,pred_df4$pred)


#5th fold

lgr5<-glm(label ~ ., family = binomial("logit"), data = f5.train)

lgr5$score<-predict(lgr5, f5.test[,2:ncol(f5.test)],type="response")
pred_df5 <- data.frame(dep_var = f5.test$label, score = lgr5$score)
pred_df5$pred<-pred_df5$score>0.5
lgr.table5<-table(pred_df5$dep_var,pred_df5$pred)


lr.perf<-ROCR::performance(prediction(lgr$score,f.test$label),'tpr','fpr')