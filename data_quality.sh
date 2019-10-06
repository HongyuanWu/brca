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
