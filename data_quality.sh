library("BART")
library("Rtsne")  
library("ggplot2")  
library("readxl")
library("BayesTree")
library("randomForest")
library("tidyverse")
library("caret")
library("MASS")


# merge bismark cov to matrix 
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/cov2matrix.pl
perl cov2matrix.pl > brca.txt

# check missing value for each cpg and each sample 
data<-read.table("brca.txt",head=T,check.names=F)
nna<-apply(data,1,function(x) sum(is.na(x)))/ncol(data)
input<-data[-which(nna>0.5),]
nna<-apply(input,2,function(x) sum(is.na(x)))/nrow(input)
sort(nna)
names(nna[nna>0.3])
write.table(names(nna[nna>0.3]),file="high.missing.sample.txt",sep="\t",quote=F)
input$BN190199
input$BN190073

data<-read.table("brca.txt",head=T,check.names=F)
nna<-apply(data,1,function(x) sum(is.na(x)))/ncol(data)
data<-data[-which(nna>0.3),]
CHR<-unlist(lapply(strsplit(rownames(data),":"),function(x) x[1]))
POS<-as.numeric(unlist(lapply(strsplit(rownames(data),":"),function(x) x[2])))
START<-POS-150
END<-POS+150
input<-data.frame(CHR,START,END)
write.table(input,file="brca.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)

cd ~/hpc/methylation/brca/19B0731C_MethylTarget/methyfreq
wget https://raw.githubusercontent.com/Shicheng-Guo/brcameth/master/extdata/brca.tcga.target.hg19.bed

bedtools sort -i brca.hg19.bed > brca.hg19.sort.bed
bedtools merge -i brca.hg19.sort.bed > brca.hg19.sort.merge.bed

bedtools intersect -wao -a brca.tcga.target.hg19.bed -b brca.hg19.sort.merge.bed | grep -v '\-1' | wc -l 
bedtools intersect -wao -a brca.tcga.target.hg19.bed -b brca.hg19.sort.merge.bed | grep '\-1' | wc -l 



# check missing value for each cpg and each sample 
data<-read.table("brca.txt",head=T,check.names=F)
nna<-apply(data,1,function(x) sum(is.na(x)))/ncol(data)
input<-data[-which(nna>0.5),]
nna<-apply(input,2,function(x) sum(is.na(x)))/nrow(input)
sort(nna)
names(nna[nna>0.3])
write.table(names(nna[nna>0.3]),file="high.missing.sample.txt",sep="\t",quote=F)
input$BN190199
input$BN190073

# differential methylation loci analysis
data<-read.table("brca.txt",head=T,check.names=F)
nna<-apply(data,1,function(x) sum(is.na(x)))/ncol(data)
input<-data[-which(nna>0.5),]
nna<-apply(input,2,function(x) sum(is.na(x)))/nrow(input)
input<-input[,-which(nna>0.3)]
system("wget https://raw.githubusercontent.com/Shicheng-Guo/brcameth/master/extdata/phen.txt")
phen<-read.table("phen.txt",head=T)

input<-input[,-which(is.na(match(colnames(input),phen[,1])))]
y<-phen[match(colnames(input),phen[,1]),3]
levels(y)<-c(0,1)
input<-data.frame(y,t(input))

Pglm<-c()
Pt<-c()
for(i in 2:ncol(input)){
temp<-na.omit(data.frame(y=input$y,x=input[,i]))
fit1 <- (bayesglm(y ~ .,family=binomial,data=temp,na.action=na.omit))
fit2<-t.test(x~y,data=temp)
Pglm<-rbind(Pglm,summary(fit)$coefficients[2,])
Pt<-rbind(Pt,c(Z=fit2$statistic,fit2$estimate,P=fit2$p.value,CI=fit2$conf.int))

CHR<-unlist(lapply(strsplit(colnames(input)[2:ncol(input)],"[.]"),function(x) x[1]))
START<-as.numeric(unlist(lapply(strsplit(colnames(input)[2:ncol(input)],"[.]"),function(x) x[2])))-150
END<-as.numeric(unlist(lapply(strsplit(colnames(input)[2:ncol(input)],"[.]"),function(x) x[2])))+150
rownames(Pt)<-colnames(input)[2:ncol(input)]
rownames(Pglm)<-colnames(input)[2:ncol(input)]
Ptt<-data.frame(CHR,START,END,Pt)
Pglmt<-data.frame(CHR,START,END,Pglm)
}
write.table(Ptt,file="Ptt.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)
write.table(Pglmt,file="Pglmt.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)


input<-data
cv.error <- NULL
k <- 5
rlt1<-c()
rlt2<-c()
V<-c()
for(i in 1:k){
  index <- sample(1:nrow(input),round(0.9*nrow(input)))
  train.cv <- input[index,]
  test.cv <- input[-index,]
  fit <- (bayesglm(Phen ~ LenMed+MethylHBV+Sex+Age,family=binomial,data=train.cv,na.action=na.omit))
  pscores <- predict(fit,test.cv)
  V=rbind(V,data.frame(test.cv$phen,pscores))
}
V
plotROC(data=V,cOutcome=1,predrisk=V$pscores)
head(V)




