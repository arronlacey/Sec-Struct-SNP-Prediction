#Protr Analysis
library(protr)
library(data.table)


# define functions

extractTCraw = function(y) {

k<-as.numeric(y[,2])
prot<-as.character(y[,1])
pos<-as.character(y[,2])
wild<-as.character(y[,3])
sub<-as.character(y[,4])

a1<-substr(y[,5],k-2,k)
a2<-substr(y[,5],k-1,k+1)
a3<-substr(y[,5],k,k+2)

xSplitted = strsplit(y[,5], split = '')[[1]]
n  = nchar(y[,5])
TC1 = summary(factor(
  paste(paste(xSplitted[-c(n, n-1)], xSplitted[-c(1, n)], sep = ''),
        xSplitted[-c(1, 2)], sep = ''),
  levels = a1), maxsum = 8001)/(n - 2)
TC2 = summary(factor(
  paste(paste(xSplitted[-c(n, n-1)], xSplitted[-c(1, n)], sep = ''),
        xSplitted[-c(1, 2)], sep = ''),
  levels = a2), maxsum = 8001)/(n - 2)
TC3 = summary(factor(
  paste(paste(xSplitted[-c(n, n-1)], xSplitted[-c(1, n)], sep = ''),
        xSplitted[-c(1, 2)], sep = ''),
  levels = a3), maxsum = 8001)/(n - 2)

raw1<-c(prot,pos,wild,sub)
raw2<-c(TC1,TC2,TC3)
rawTC<-c(raw1,raw2)
rawnum<-as.numeric(raw2[-grep("'", names(raw2))])
c(raw1,rawnum)

}


extractCTDCraw = function(x) {
  
  k<-as.numeric(x[,2])
  prot<-as.character(x[,1])
  pos<-as.numeric(x[,2])
  wild<-as.character(x[,3])
  sub<-as.character(x[,4])
  
  group1 = list(
    'hydrophobicity'  = c('R', 'K', 'E', 'D', 'Q', 'N'),
    'normwaalsvolume' = c('G', 'A', 'S', 'T', 'P', 'D', 'C'),
    'polarity'        = c('L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'),
    'polarizability'  = c('G', 'A', 'S', 'D', 'T'),
    'charge'          = c('K', 'R'),
    'secondarystruct' = c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'),
    'solventaccess'   = c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'))
  
  group2 = list(
    'hydrophobicity'  = c('G', 'A', 'S', 'T', 'P', 'H', 'Y'),
    'normwaalsvolume' = c('N', 'V', 'E', 'Q', 'I', 'L'),
    'polarity'        = c('P', 'A', 'T', 'G', 'S'),
    'polarizability'  = c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'),
    'charge'          = c('A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L',
                          'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
    'secondarystruct' = c('V', 'I', 'Y', 'C', 'W', 'F', 'T'),
    'solventaccess'   = c('R', 'K', 'Q', 'E', 'N', 'D'))
  
  group3 = list(
    'hydrophobicity'  = c('C', 'L', 'V', 'I', 'M', 'F', 'W'),
    'normwaalsvolume' = c('M', 'H', 'K', 'F', 'R', 'Y', 'W'),
    'polarity'        = c('H', 'Q', 'R', 'K', 'N', 'E', 'D'),
    'polarizability'  = c('K', 'M', 'H', 'F', 'R', 'Y', 'W'),
    'charge'          = c('D', 'E'),
    'secondarystruct' = c('G', 'N', 'P', 'S', 'D'),
    'solventaccess'   = c('M', 'S', 'P', 'T', 'H', 'Y'))
  
  xSplitted = substr(x[1,5],pos,pos)
  
  # Get groups for each property & each amino acid
  
  g1 = lapply(group1, function(g) which(xSplitted %in% g))
  names(g1) = paste(names(g1), '1.', sep = '.')
  g2 = lapply(group2, function(g) which(xSplitted %in% g))
  names(g2) = paste(names(g2), '2.', sep = '.')
  g3 = lapply(group3, function(g) which(xSplitted %in% g))
  names(g3) = paste(names(g3), '3.', sep = '.')

  CTDC1 <- unlist(c(g1, g2, g3))
  CTDC1.names <- names(CTDC1)
  CTDC2 <- as.numeric(sub("\\D*(\\d+).*", "\\1",names(CTDC1[order(as.character(CTDC1.names))])))
  CTDC <- c(prot,pos,wild,sub,CTDC2)
  CTDC
}

####################################################################################




#load in dataset
t<-read.csv("~/Phd/scripts/humvar3/humvar4protr.csv", stringsAsFactors=FALSE)
names(t)<-c("prot","pos","wild","sub","seq")

# alternate sequences with wild subbed in

w<-t
# sub in wild type to sequence
w[,5]<-paste0(substr(w[,5],1,as.numeric(w[,2])-1),w[,3],substr(w[,5],as.numeric(w[,2])+1,nchar(w[,5])))


####### extract tripeptide composition frequencies found around the SNP

TC<-lapply(1:nrow(t), function(x) extractTCraw(t[x,1:5]))   # specify window of 2 amino acids either side of mutation
TC.wild<-lapply(1:nrow(w), function(x) extractTCraw(w[x,1:5]))


####### Extract physiochemical properties at SNP position

CTDC<-lapply(1:nrow(t), function(x) extractCTDCraw(t[x,1:5]))
CTDC.wild<-lapply(1:nrow(w), function(x) extractCTDCraw(w[x,1:5]))

# merge physiochemical properties with tripeptide frequencies  

protr.mut<-mapply(c, TC, CTDC, SIMPLIFY = F)
protr.mut.clean<-protr.mut[sapply(protr.mut, length) == 18]
protr.wild<-mapply(c, TC.wild, CTDC.wild, SIMPLIFY = F)
protr.wild.clean<-protr.wild[sapply(protr.wild, length) == 18]

testmut<-do.call(rbind.data.frame, protr.mut)
testmut.clean<-do.call(rbind.data.frame, protr.mut.clean)
names(testmut)<-c("prot","pos","wild","mut","freq.win1","freq.win2","freq.win3","prot2","pos2","wild2","mut2","charge","hydrophobicity","normwaalsvolume","polarity","polarizability","secondarystruct","solventaccess")
names(testmut.clean)<-c("prot","pos","wild","mut","freq.win1","freq.win2","freq.win3","prot2","pos2","wild2","mut2","charge","hydrophobicity","normwaalsvolume","polarity","polarizability","secondarystruct","solventaccess")
testmut<-subset(testmut, select = -c(prot2,pos2,wild2,mut2))
testmut.clean<-subset(testmut.clean, select = -c(prot2,pos2,wild2,mut2))
cols = c(2, 5, 6, 7);    
testmut.clean[,cols] = apply(testmut.clean[,cols], 2, function(x) toint(x));

testwild<-do.call(rbind.data.frame, protr.wild)
testwild.clean<-do.call(rbind.data.frame, protr.wild.clean)
names(testwild)<-c("prot","pos","wild","mut","freq.win1","freq.win2","freq.win3","prot2","pos2","wild2","mut2","charge","hydrophobicity","normwaalsvolume","polarity","polarizability","secondarystruct","solventaccess")
names(testwild.clean)<-c("prot","pos","wild","mut","freq.win1","freq.win2","freq.win3","prot2","pos2","wild2","mut2","charge","hydrophobicity","normwaalsvolume","polarity","polarizability","secondarystruct","solventaccess")
testwild<-subset(testwild, select = -c(prot2,pos2,wild2,mut2))
testwild.clean<-subset(testwild.clean, select = -c(prot2,pos2,wild2,mut2))
cols = c(2, 5, 6, 7);    
testwild.clean[,cols] = apply(testwild.clean[,cols], 2, function(x) toint(x));


protr_clean<-unique(merge(testmut.clean,testwild.clean,by.x=c("prot","pos","wild","mut"),
                   by.y=c("prot","pos","wild","mut"),all.x=TRUE,all.y=TRUE))
protr<-protr_clean

#protr<-unique(merge(mut=testmut,wild=testwild,by.x=c("prot","pos","wild","mut"),
 #                       by.y=c("prot","pos","wild","mut"),all.x=TRUE,all.y=TRUE,stringsAsFactors = FALSE))

# replaces . in column with _
#names(protr)<-gsub(x = names(protr), pattern = "\\.", replacement = "_")
#names(protr.clean)<-gsub(x = names(protr_clean), pattern = "\\.", replacement = "_")


#protr[,c(5:7,19:21)] <- sapply(protr[,c(5:7,19:21)],toint)
#protr<-protr[- grep("script", protr),]
#protr<-factor(protr)


