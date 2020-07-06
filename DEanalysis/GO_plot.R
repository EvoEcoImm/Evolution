#! /usr/bin/Rscript
#cockroach#
library(readr)
gos<- read.csv("tgoplot",header=F,sep="\t")
gos$V2<- -log10(gos$V8)
gos$V3<- gos$V4/gos$V5
gobps<- gos[gos$V7=="BP",]
gomfs<- gos[gos$V7=="MF",]
'''dotplot'''
plot(gomfs$V2 ~ gomfs$V3, cex=log2(gomfs$V4),xlab="Proportation of regulated genes",ylab="-log10(adjust p-value)",bty="l",cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="blue",las=1,xlim=c(0,0.7),ylim=c(1.3,2.3))
points(gobps$V2 ~ gobps$V3, cex=log2(gobps$V4),cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="red")
legend(0.60,2.0, legend=c("BP","MF"), col=c("red","blue"),inset=c(-0.2,0.2),pch=15, pt.cex=4,title="GO Category",bty="n",y.intersp=1.5,x.intersp=2,cex=1.5,xpd=T)
legend(0.5,1.7, legend=c("5","10","20"), inset=c(-0.2,0.2),pch=1, pt.cex=c(2.3,3.3,4.3),title="Number of regulated genes",bty="n", y.intersp= 1.5,x.intersp=1.5, cex=1.5,xpd=T)
dev.print(pdf,"tGO_plot.pdf")
write_csv(gos,"./tGO_plot.csv")

library(readr)
gos<- read.csv("cgoplot",header=F,sep="\t")
gos$V2<- -log10(gos$V8)
gos$V3<- gos$V4/gos$V5
gobps<- gos[gos$V7=="BP",]
gomfs<- gos[gos$V7=="MF",]
'''dotplot'''
plot(gomfs$V2 ~ gomfs$V3, cex=log2(gomfs$V4),xlab="Proportation of regulated genes",ylab="-log10(adjust p-value)",bty="l",cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="blue",las=1,xlim=c(0,0.4),ylim=c(1,15))
points(gobps$V2 ~ gobps$V3, cex=log2(gobps$V4),cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="red")
legend(0.3,14, legend=c("BP","MF"), col=c("red","blue"),inset=c(-0.2,0.2),pch=15, pt.cex=4,title="GO Category",bty="n",y.intersp=1.5,x.intersp=2,cex=1.5,xpd=T)
legend(0.25,10, legend=c("5","10","20"), inset=c(-0.2,0.2),pch=1, pt.cex=c(2.3,3.3,4.3),title="Number of regulated genes",bty="n", y.intersp= 1.5,x.intersp=1.5, cex=1.5,xpd=T)
dev.print(pdf,"cGO_plot.pdf")
write_csv(gos,"./cGO_plot.csv")


q()

#termite#
#reproductives
library(readr)
rtgos<- read.csv("rtreat.goseq_BP.cut",header=F,sep="\t")
rtgos$V2<- -log10(rtgos$V8)
rtgos$V3<- rtgos$V4/rtgos$V5
plot(rtgos$V2 ~ rtgos$V3, cex=log2(rtgos$V4),xlab="Proportation of regulated genes",ylab="-log10(adjust p-value)",bty="l",cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="red",las=1,pch=19,xlim=c(0,0.4),ylim=c(1.2,5))
#soldiers
gos<- read.csv("tsgoplot",header=F,sep="\t")
gos$V2<- -log10(gos$V8)
gos$V3<- gos$V4/gos$V5
gobps<- gos[gos$V7=="BP",]
gomfs<- gos[gos$V7=="MF",]
points(gobps$V2 ~ gobps$V3, cex=log2(gobps$V4),cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="red")
points(gomfs$V2 ~ gomfs$V3, cex=log2(gomfs$V4),cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="blue")
#workers
wtgos<- read.csv("wtreat.goseq_BP.cut",header=F,sep="\t")
wtgos$V2<- -log10(rtgos$V8)
wtgos$V3<- wtgos$V4/wtgos$V5
points(wtgos$V2 ~ wtgos$V3, cex=log2(wtgos$V4),cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="red",pch=10)
wtgos<- read.csv("wtreat.goseq_MF.cut",header=F,sep="\t")
wtgos$V2<- -log10(rtgos$V8)
wtgos$V3<- wtgos$V4/wtgos$V5
points(wtgos$V2 ~ wtgos$V3, cex=log2(wtgos$V4),cex.lab=2,cex.axis=1.5,tcl=0.3,lwd=2, col="blue",pch=10)


legend(0.33,5.0, legend=c("BP","MF"), col=c("red","blue"),inset=c(-0.2,0.2),pch=15, pt.cex=3,title="GO Category",bty="n",y.intersp=1,x.intersp=2,cex=1.5,xpd=T)
legend(0.3,4, legend=c("5","10","20"), inset=c(-0.2,0.2),pch=1, pt.cex=c(2.3,3.3,4.3),title="Number of regulated genes",bty="n", y.intersp= 1.5,x.intersp=1.5, cex=1.5,xpd=T)
legend(0.2,5,legend=c("Reproductives","Soldier","Worker"),pch=c(19,1,10), bty="n",inset=c(-0.2,0.2), pt.cex=3,title="Caste",y.intersp= 1,x.intersp=1.5, cex=1.5,xpd=T)
dev.print(pdf,"terGO_plot.pdf")
	

