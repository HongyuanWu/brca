
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

