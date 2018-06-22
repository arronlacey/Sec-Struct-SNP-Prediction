toint <- function(x){as.numeric(as.character(x))}

library(sqldf)
library(RSQLite)
library(RWeka)
library(foreign)
library(randomForest)
library(stringr)
library(rpart)
library(caret)
library(SuperLearner)
library(plyr)
library(e1071)
library(stats)

# set up sqlite database
db <- dbConnect(RSQLite::SQLite(), 'Test1.sqlite')

# load base data

humvar <- read.csv('/Users/arron/Downloads/humvar3.csv',sep=',',row.names=NULL, stringsAsFactors = FALSE, na.strings = "NA")
humvar_psipred <- read.csv('/Users/arron/Documents/Phd/scripts/humvar3/merged.csv',sep=',',row.names=NULL, stringsAsFactors = FALSE, na.strings = "NA")
vep <- read.csv('/Users/arron/Documents/Phd/scripts/humvar3/vep_38.txt',sep='\t',stringsAsFactors = FALSE, na.strings = "NA")
source("Phd/scripts/protr/protr.r")

# rename some columns

vep$Amino_acids<-str_split_fixed(vep$Amino_acids, "/", 2)
vep$aawild<-vep$Amino_acids[,1]
vep$aamut<-vep$Amino_acids[,2]
vep<-vep[-c(18)]
vep<-vep[!is.na(vep$PolyPhen),]
vep<-vep[!is.na(vep$SIFT),]
vep$PolyPhen_score<-as.numeric(unlist(regmatches(vep$PolyPhen,
                                                 gregexpr("[[:digit:]]+\\.*[[:digit:]]*",vep$PolyPhen))
)      )
vep$SIFT_score<-as.numeric(unlist(regmatches(vep$SIFT,
                                             gregexpr("[[:digit:]]+\\.*[[:digit:]]*",vep$SIFT))
)      )
vep$PolyPhen_pred<-regmatches(vep$PolyPhen, regexpr("^[^(]*", vep$PolyPhen))
vep$SIFT<-regmatches(vep$SIFT, regexpr("^[^(]*", vep$SIFT))
vep$PolyPhen<-NULL
vep$SIFT<-NULL

# store in database

dbWriteTable(db, 'humvar', humvar, row.names = FALSE,overwrite=TRUE)
dbWriteTable(db, 'humvar_psipred', humvar_psipred, row.names = FALSE,overwrite=TRUE)
dbWriteTable(db, 'vep', vep, row.names = FALSE,overwrite=TRUE)
dbWriteTable(db, 'protr', protr, row.names = FALSE,overwrite=TRUE)

# join the 4 tables

merged<-sqldf("select vep.*, humvar.*, protr.*,humvar_psipred.diff1,humvar_psipred.diff2,humvar_psipred.diff3, ss_x as mut_sec, ss_y as wild_sec
      from vep inner join humvar
              on vep.SWISSPROT = humvar.Prot
              and vep.Protein_position = humvar.X
              and vep.aawild = humvar.aawild
              and vep.aamut = humvar.aamut
              join humvar_psipred
              on humvar.Prot = humvar_psipred.prot
              and humvar.X = humvar_psipred.pos
              and humvar.aawild = humvar_psipred.wild
              and humvar.aamut = humvar_psipred.mut
              inner join protr
      on protr.prot = vep.SWISSPROT
              and vep.Protein_position = protr.pos
              and vep.aawild = protr.wild
              and vep.aamut = protr.mut")

# Remove "unclassified SNPs"

merged<-subset(merged, label != "Unclassified")
merged<-subset(merged, label != "Unclassifie")


# clean merged data ready for machine learning

f<-merged

f[f == "-"] <- NA    # convert dashes to NA 
f<-f[!grepl("Unclassified", f$label),]
f$label <- as.factor(f$label)

# change NA freqencies to 0.000001
f.freq<-f[,grepl("AF$", names(f))] # create dataframe of frequencies for f
f<-f[,grep("AF$", names(f),invert=TRUE)]  # remove frequencies from f
f.freq[is.na(f.freq)] <- 0.000001
f.freq<-lapply(f.freq, function(x) gsub("NA", "-", x))
f.freq <- as.data.frame(lapply(f.freq, function(x) toint(x)))
freq<-cbind(f,f.freq)

df<-freq
# check difference between pyhsiochemical states between wild and mutation
df$hydrophobicity<-df$hydrophobicity.x == df$hydrophobicity.y
df$charge<-df$charge.x == df$charge.y
df$normwaalsvolume<-df$normwaalsvolume.x == df$normwaalsvolume.y
# removes redundant columns from merge
df<-df[, !grepl("\\.x", colnames(df))]
df<-df[, !grepl("\\.y", colnames(df))]
# keep scores as features
df.score<-df[, grep("_score", colnames(df))]
df.pred<-df[, grep("_pred|phred", colnames(df))]
# these scores I'm not score how to interpret, exclude.
df.score$MutPred_score<-NULL
df.score$MNACAP_score<-NULL
df.score$LRT_score<-NULL
# keep rankscores as features
df.rankscore<-df[, grep("_rankscore", colnames(df))]
df.rankscore$MutPred_rankscore<-NULL
df.rankscore$MNACAP_rankscore<-NULL
df.rankscore$LRT_converted_rankscore<-NULL
# some misc columns kept as features
df.misc<-df[,c("Reliability_index","GERP.._NR","SiPhy_29way_logOdds","SiPhy_29way_logOdds_rankscore",
               "phastCons20way_mammalian","charge","hydrophobicity","normwaalsvolume","diff1","diff2","diff3","wild_sec","mut_sec")]
# frequency as features i.e. gnomad
df.freq<-df[,grepl("AF$", names(df))]  
df2<-cbind.data.frame(df$label,df.score,df.rankscore,df.misc,df.freq$gnomAD_AF,df.freq$gnomAD_exomes_AF,df.freq$gnomAD_genomes_AF) #orig column selection
# secondary strucutre features
df2$mut_sec<-as.factor(df2$mut_sec)
df2$wild_sec<-as.factor(df2$wild_sec)
colnames(df2)[1]<-"label"
# remove rows containing NULLs
df2<-df2[complete.cases(df2),]
df2<-df2[, !sapply(df2, is.character)]
df2 <- df2[, !duplicated(colnames(df2))]
# add gnomad
colnames(df2)[52]<-"gnomAD_AF"
colnames(df2)[53]<-"gnomAD_exomes_AF"
colnames(df2)[54]<-"gnomAD_genomes_AF"

source("Phd/Thesis/code/rfcrossfold.r")
source("Phd/Thesis/code/nbcrossfold.r")
source("Phd/Thesis/code/nncrossfold.r")
source("Phd/Thesis/code/lgrcrossfold.r")
source("Phd/Thesis/code/svmcrossfold.r")
source("Phd/Thesis/code/dtcrossfold.r")


ConfusionResults<-function(t){
  sen<-t[1]/(t[1]+t[3])
  spec<-t[4]/(t[4]+t[2])
  acc<-(t[1]+t[4])/(t[1]+t[2]+t[3]+t[4])
  tp<-t[1]
  fp<-t[3]
  fn<-t[2]
  tn<-t[4]
  disease<-t[1]+t[2]
  poly<-t[3]+t[4]
  tline<-paste0(disease,",",poly,",",tp,",",fp,",",fn,",",tn,",",
                100*round(sen,4),",",100*round(spec,4),",",100*round(acc,4))
  return(tline)
}


confusiontabs<-list(rf.table1,rf.table2,rf.table3,rf.table4,rf.table5,
                    lgr.table1,lgr.table2,lgr.table3,lgr.table4,lgr.table5,
                    nn.table1,nn.table2,nn.table3,nn.table4,nn.table5,
                    nb.table1,nb.table2,nb.table3,nb.table4,nb.table5,
                    svm.table1,svm.table2,svm.table3,svm.table4,svm.table5,
                    dt.table1,dt.table2,dt.table3,dt.table4,dt.table5
)

crossfold.table<-unlist(lapply(1:length(confusiontabs), 
                               function(i) ConfusionResults(unlist(confusiontabs[i]))))
read.table(text = crossfold.table, sep = rawToChar(as.raw(2)), header = FALSE
)




f2.test.4roc<-f4.test[ , grepl( "_rankscore" , names( f4.test ) ) ]    # toggle between rankscore and score
f2.test.4roc<-f2.test.4roc[ , !grepl( "_fitCons" , names( f2.test.4roc ) ) ]
f2.test.4roc<-f2.test.4roc[ , !grepl( "converted" , names( f2.test.4roc ) ) ]
f2.test.4roc$Polyphen<-f4.test$PolyPhen_score
f2.test.4roc<-cbind(f4.test$label,f2.test.4roc)
f2.test.4roc<-f2.test.4roc[, -grep("1$", colnames(f2.test.4roc))]

#h.pd<-ROCR::prediction(f2.test$CADD_raw_rankscore, f2.test$label)
#table(observed = f2.test$label, predicted = h.pd)

#hum.perf<-ROCR::performance(ROCR::prediction(1-f2.test.4roc$REVEL_rankscore,f2.test$label),'tpr','fpr',measure="auc")
f2.test.4roc$SIFT<-1-f4.test$SIFT_score


hum.perfs <- lapply(2:ncol(f2.test.4roc), function(i) 
  ROCR::performance(ROCR::prediction(1-f2.test.4roc[,i],f2.test.4roc[,1]),'tpr','fpr'))

## AUC

hum.perfs.auc <- lapply(2:ncol(f2.test.4roc), function(i) 
  ROCR::performance(ROCR::prediction(1-f2.test.4roc[,i],f2.test.4roc[,1]),'tpr','fpr',measure="auc"))

auc<-unlist(lapply(seq_along(hum.perfs.auc), function(x) hum.perfs.auc[[x]]@y.values))
# attach classifiers to AUC value and sort
auc2<-cbind.data.frame(auc,names(f2.test.4roc[,2:length(f2.test.4roc)]))
auc3<-auc2[order(-auc),]
# extract classifier names
auc3<-auc3$`names(f2.test.4roc[, 2:length(f2.test.4roc)])`
auc3<-as.character(auc3)
#clean up classifier names
auc.names<-gsub("\\_.*","",auc3)

hum.perfs.auc
# set roc curve color in order of classifier plotted
roc.col<-distinctColorPalette(nrow(auc2))
col2<-cbind.data.frame(auc2,roc.col)
col3<-col2[order(-auc),]
col3<-as.character(col3$roc.col)


# plot ROC curves

# humvar_comparison.png

# space for legend
par(mar=c(5.1, 4.1, 4.1, 10.5), xpd=TRUE)
# plot my classifier first
plot(rf.perf, col='black',lwd="1.2",main="Humvar Test Set Classification") # my classifier
# add further classifier plots
for (i in 1:length(hum.perfs)){
  plot(hum.perfs[[i]], col=roc.col[i],add=TRUE,lwd="1.2") 
}
# legend
hum.perf.names<-gsub("\\_.*","",names(f2.test.4roc[,2:ncol(f2.test.4roc)]))
legend("topright", inset=c(-0.40,0),
       c("Random Forest",auc.names),
       col=c("black",col3),
       lty=c(1,1),
       bty="n",
       cex=0.8)





############################################# get cut off for 0.9 fpr

getROCStats<-function(k,i){
  r<-ROCR::performance(ROCR::prediction(1-k[,i],k[,1]),'tpr','fpr')
  r.auc<-ROCR::performance(ROCR::prediction(1-k[,i],k[,1])
                           ,'tpr','fpr',measure="auc")
  cutoffs <- data.frame(cut=r@alpha.values[[1]], fpr=r@x.values[[1]], 
                        tpr=r@y.values[[1]])
  cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
  stats<-cbind(names(k[i]),head(subset(cutoffs, tpr < 0.95),n=1),unlist(r.auc@y.values))
  names(stats)<-c("Classifier","Threshold","FPR","TPR","AUC")
  stats$Threshold<-round(stats$Threshold,2)
  stats$FPR<-round(stats$FPR,2)
  stats$AUC<-round(stats$AUC,2)
  stats
}             


getROCStats(rf4,2)

rocstat<-c()
for (i in 2:ncol(f2.test.4roc)){
  rocstat[[i]]<-getROCStats(f2.test.4roc,i)
}

roc.table<-ldply(rocstat, data.frame)
roc.table<-roc.table[order(-roc.table$AUC),]
roc.table<-roc.table[,-2]
roc.table$Classifier<-gsub("\\_.*","",roc.table$Classifier)
roc.tab<-rbind(roc.table,c("Random Forest",stats$FPR,stats$TPR,stats$AUC))
roc.tab<-roc.tab[order(roc.tab$FPR),]
roc.tab<-roc.tab[roc.tab$TPR >=0.93,]
roc.tab<-data.frame(roc.tab$Classifier,roc.tab$FPR)
roc.tab<-roc.tab[!duplicated(roc.tab), ]
names(roc.tab)<-c("Classifier","FPR")
roc.tab$Classifier <- factor(roc.tab$Classifier, levels = roc.tab$Classifier[order(roc.tab$FPR)])

roc.tab$FPR<-toint(roc.tab$FPR)

ggplot(roc.tab,aes(x=Classifier,y=FPR)) + geom_bar(stat="identity") +
  theme_bw() +
  ggtitle("Specificty plot at 95% sensitivity") +
  xlab("Classifier") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous()


