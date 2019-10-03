# merge bismark cov to matrix 
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/cov2matrix.pl
perl cov2matrix.pl > brca.txt

# check missing value for each cpg and each sample 
data<-read.table("brca.txt",head=T,check.names=F)
nna<-apply(data,1,function(x) sum(is.na(x)))/ncol(data)
input<-data[-which(nna>0.5),]
nna<-apply(input,2,function(x) sum(is.na(x)))/nrow(input)
sort(nna)
input$BN190199
input$BN190073
