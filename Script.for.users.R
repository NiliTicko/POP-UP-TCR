
#R Packages to use:

library('reshape2')
library('stringr')
library('randomForest')

#Your dataset file (named here 'tcr_dataset') should have the following columns:
#A column of CDR3 beta sequences, termed here column 'beta' (tcr_dataset$beta)
#A column of peptide sequences, termed here 'ag'
#A column of TCR identification number, termed here 'tcr'
#A column of status, e.g., 'healthy/ sick', termed here 'status'


#Read your dataset

tcr_dataset= read.csv("*path/file.to.read.csv")

#Make sliding windows from the peptide sequences:

tcr_data= tcr_dataset

tcr_data$ag= paste("XX",tcr_data$ag,sep="")
tcr_data$ag= paste(tcr_data$ag,"XX",sep="")

	for(i in(1:length(unique(tcr_data$tcr)))){
a= subset(tcr_data, tcr_data$tcr== unique(tcr_data$tcr)[i])
a$letters=str_split_fixed(a$ag, "", n=nchar(a$ag))

r= data.frame(a$letters)
	for(p in(1:(length(r)-4))){
d3=data.frame(unique(a$tcr), r[ ,(p):(p+4)],p, unique(a$status))
#write.table(d3, file = "*path/file_ag.csv", sep = ",", append = TRUE, quote = FALSE,  col.names = FALSE, row.names = FALSE)
write.table(d3, file = "C:/Nili/Pred_clust/test/file_ag.csv", sep = ",", append = TRUE, quote = FALSE,  col.names = FALSE, row.names = FALSE)

}
}

ag.wind= read.csv("*path/file_ag.csv", header=F)

colnames(ag.wind)= c("tcr", "ag1","ag2","ag3","ag4","ag5","loc.ag","status")

#Make sliding windows from the CDR3 beta sequences: 

tcr_data$b= paste("XX",tcr_data$beta,sep="")
tcr_data$b= paste(tcr_data$b,"XX",sep="")

	for(i in(1:length(unique(tcr_data$tcr)))){
a= subset(tcr_data, tcr_data$tcr== unique(tcr_data$tcr)[i])
a$letters=str_split_fixed(a$b, "", n=nchar(a$b))

r= data.frame(a$letters)
	for(p in (1:(length(r)-4))){
d1=data.frame(unique(a$tcr), r[,(p):(p+4)],p)
#write.table(d1, file = file_beta, sep = ",", append = TRUE, quote = FALSE,
#  col.names = FALSE, row.names = FALSE)
write.table(d1, file = "*path/beta.csv", sep = ",", append = TRUE, quote = FALSE,
  col.names = FALSE, row.names = FALSE)

}
}
tb.wind= read.csv("*path/beta.csv", header=F)
colnames(tb.wind)= c("tcr", "aa1","aa2","aa3","aa4","aa5","loc.aa")

#Prepare the vectors for the 1st machine:

aa=as.array(c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X'))

a.cod=matrix(0, nrow=21, ncol=21)
rownames(a.cod)=unique(aa)
colnames(a.cod)=1:21
for (i in (1: 21)){
a.cod[i,i]=1
}

aa.code=as.data.frame(a.cod)
aa.code$AA=rownames(aa.code)
aa.code$Cod=paste(c(aa.code[,1:21]))
aa.code$code= gsub(",","", aa.code$Cod)

aa.code$code= gsub("[c][(]","", aa.code$code)
aa.code$code= gsub("[)]","", aa.code$code)


an=as.array(1:21)
an.cod=matrix(0, nrow=21, ncol=21)
rownames(an.cod)=unique(an)
colnames(an.cod)=1:21
for (i in (1: 21)){
an.cod[i,i]=1
}

an.code=as.data.frame(an.cod)
an.code$Cod=paste(c(an.code[,1:21]))
an.code$code= gsub(",","", an.code$Cod)

an.code$code= gsub("[c][(]","", an.code$code)
an.code$code= gsub("[)]","", an.code$code)
an.code$Num=rownames(an.code)
an.code$Cod=paste(c(an.code[,1:21]))

merged= merge(tb.wind, ag.wind, by='tcr') 
data1= merged[,c(1:6,8:12,7,13,14)]

data1$num.aa= data1$loc.aa
data1$num.ag= data1$loc.ag

for (i in (2:(ncol(data1)-5))){
ind=match(data1[,i], aa.code$AA)
data1[,i]=aa.code$code[ind]
}

ind1= match(data1$loc.aa, an.code$Num)
data1$loc.aa= an.code$code[ind1]
ind2= match(data1$loc.ag, an.code$Num)
data1$loc.ag= an.code$code[ind2]

for(n in(2:11)){
data1[,n]= factor(data1[,n], levels= levels(factor(aa.code$code)))
}
 data1$loc.aa= factor(data1$loc.aa, levels= levels(factor(an.code$code)))
 data1$loc.ag= factor(data1$loc.ag, levels= levels(factor(an.code$code)))

#Read the model of the 1st machine(Attached in the file) :

RF1= readRDS("*path/rf.beta.1st.rds")

beta.pred= predict(RF1, data1[,2:13], type="prob")

pred.b= cbind(data1, beta.pred)

#Prepare data predicted by 1st machine for the 2nd machine:

pred.b$score= round(pred.b[,18], digits=1)
pred2= pred.b[, c(1,14,15,16,19)]

pred2= subset(pred2, pred2$num.ag < 16)
pred2= subset(pred2, pred2$num.aa < 14)

#Take the scores for each TCR: order by peptide locations

	for(m in(1:length(unique(pred2$tcr)))){
tcr1= subset(pred2, pred2$tcr== unique(pred2$tcr)[m])

#Sort by ag location and also by cdr aa location:

tcr2= tcr1[order(tcr1$num.ag, tcr1$num.aa),]

#Make a matrix of ag. location by cdr aa:

mat.bet= matrix(0.00001, nrow= 13, ncol= 15)
rownames(mat.bet)= paste(1:13,"cdr", sep="")
colnames(mat.bet)= paste(1:15, "ag", sep="")

#Put all scores from TCR by ag. location per cdr aa

	for(i in(1:length(unique(tcr2$num.ag)))){
ag= subset(tcr2, tcr2$num.ag == unique(tcr2$num.ag)[i])
	for(n in(1:length(unique(ag$num.aa)))){
aa= subset(ag, ag$num.aa== unique(ag$num.aa)[n])
mat.bet[n,i]= aa$score

dmat= as.data.frame(mat.bet)

#Line up all matrix contents into vectors for the 2nd machine:

score.in.line= c(dmat[,1], dmat[,2], dmat[,3], dmat[,4],dmat[,5],dmat[,6], dmat[,7], dmat[,8], dmat[,9],dmat[,10],dmat[,11], dmat[,12], dmat[,13], dmat[,14],dmat[,15])
}
}
write.table(data.frame(unique(aa$tcr), unique(aa$status), as.list(score.in.line)), file = "*path/vector.for.2nd.csv", sep = ",", append = TRUE, quote = FALSE,
col.names = FALSE, row.names = FALSE)
}

test.set= read.csv("*path/vector.for.2nd.csv", header= F)

b195= paste('b', 1:195, sep='')
colnames(test.set)= c('tcr', 'status', b195)

#Read the model for the 2nd machine (Attached in the file):

RF2= readRDS("*path/rf.only.b.2nd.rds")

pred= predict(RF2, test.set[,3:197], type= "prob") 
pred_test_full= cbind(test.set, pred)

pred_test= pred_test_full[,c(1,2,199)]
colnames(pred_test)= c("tcr","status","predictionScore")
 
indi= match(pred_test$tcr, tcr_dataset$tcr)
pred_test$CDR3.beta= tcr_dataset$beta[indi]
pred_test$peptide= tcr_dataset$ag[indi]

#Write the file with the prediction for each CDR# beta-peptide pair:

write.csv(pred_test, "*path/tcr.pep.predicted.csv")