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
cd ~/hpc/methylation/brca/19B0731C_MethylTarget/methyfreq
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/cov2matrix.pl
perl cov2matrix.pl > brca.txt

######## methyPlot ################
mkdir methplot
awk -F '[\t:]' '{print $1,$2,$2+1}' OFS="\t" brca.txt | grep chr | bedtools intersect -wao -a - -b brca.tcga.target.hg19.bed | grep -v '\-1' | awk '{print $1":"$2,$7}' > brca.map

data<-read.table("brca.txt",head=T,check.names=F)
map<-read.table("brca.map",sep="")
phen<-read.table("phen.txt",head=T)

for(i in unique(map[,2])){
temp<-t(data[match(map[map[,2] %in% i,1],rownames(data)),])
class<-as.character(phen[match(rownames(temp),phen[,1]),3])
position=unlist(lapply(strsplit(colnames(temp),":"),function(x) x[2]))
temp<-cbind(temp,class)
temp<-rbind(temp,position)
temp[nrow(temp),ncol(temp)]<-""
write.table(temp,file=paste("./methplot/",i,".methploter.input.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=T)
}
######## methyPlot ################


### check missing value for each cpg and each sample 
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
Pglm<-rbind(Pglm,summary(fit1)$coefficients[2,])
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

bedtools intersect -wao -a Ptt.hg19.bed -b brca.tcga.target.hg19.bed > Ptt.intersect.hg19.bed
bedtools intersect -wao -a Pglmt.hg19.bed -b brca.tcga.target.hg19.bed > Pglmt.intersect.hg19.bed

## R plot for dual-y-axis figure
setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/brca/19B0731C_MethylTarget/methyfreq")
data<-read.table("Ptt.intersect.hg19.bed")
head(data)
for(i in unique(data$V13)){
  png(paste(i,".png",sep=""))
  temp<-subset(data,V13==i)
  par(mar = c(5,5,2,5))
  with(temp,plot(x=V3,y=-V4,type="l",main=i,xlab="Genomic Position (bp)",lwd=2,ylab="delta beta % (T-N)",col="red"))
  par(new = T)
  with(temp,plot(x=V3,y=-log(V7,10),pch=16, xlab=NA, cex=1.5,col="blue",ylab=NA,axes=F,ylim=c(0,max(-log(V7,10)))))
  axis(side = 4)
  mtext(side = 4, line = 3, expression(-log[10](italic(p))))
  legend("bottomright",legend=c( "delta beta % (T-N)",expression(-log[10](italic(p)))),cex=1.5,lwd=2,lty=c(1,0), pch=c(NA, 16), col=c("red", "blue"),bty="n")
  dev.off()
  print(paste(i,mean(-temp$V4),mean(-log(temp$V7,10)),sep=" "))
}




