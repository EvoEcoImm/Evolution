#! /usr/bin/Rscript

library(dplyr)
boup<-read.csv("/media/shulinhe/DATA/Full_analysis/Bo_full/bo_salmon/bo_ii_TvsC_sig_up.csv",header=T,sep=",")
boim<-read.csv("/media/shulinhe/DATA/Full_analysis/immune/bo_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(boim)<- c("clust","Component.ID","Family.Name","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
boimup<- boim[boim$Component.ID %in% boup$row,]
boimupcount<- as.data.frame(table(boimup$Family.Name))
colnames(boimupcount)<- c("Family","boup")
boimcount<- as.data.frame(table(boim$Family.Name))
colnames(boimcount)<- c("Family","bototal")

cmup<-read.csv("/media/shulinhe/DATA/Full_analysis/Cm_full/cm_salmon/cmres_TvsC_sig_up.csv",header=T,sep=",")
cmim<-read.csv("/media/shulinhe/DATA/Full_analysis/immune/cm_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(cmim)<- c("clust","Component.ID","Family.Name","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
cmimup<- cmim[cmim$Component.ID %in% cmup$row,]
cmimupcount<- as.data.frame(table(cmimup$Family.Name))
colnames(cmimupcount)<- c("Family","cmup")
cmimcount<-as.data.frame(table(cmim$Family.Name))
colnames(cmimcount)<- c("Family","cmtotal")

rup<-read.csv("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/rres_TvsC_sig_up.csv",header=T,sep=",")
sup<-read.csv("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/sres_TvsC_sig_up.csv",header=T,sep=",")
wup<-read.csv("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/wres_TvsC_sig_up.csv",header=T,sep=",")
ncim<-read.csv("/media/shulinhe/DATA/Full_analysis/immune/nc_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(ncim)<- c("clust","Component.ID","Family.Name","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
rimup<- ncim[ncim$Component.ID %in% rup$row,]
simup<- ncim[ncim$Component.ID %in% sup$row,]
wimup<- ncim[ncim$Component.ID %in% wup$row,]
rimupcount<- as.data.frame(table(rimup$Family.Name))
colnames(rimupcount)<- c("Family","rup")
simupcount<- as.data.frame(table(simup$Family.Name))
colnames(simupcount)<- c("Family","sup")
wimupcount<- as.data.frame(table(wimup$Family.Name))
colnames(wimupcount)<- c("Family","wup")
ncimcount<- as.data.frame(table(ncim$Family.Name))
colnames(ncimcount)<- c("Family","nctotal")
Sumtable<- full_join(boimcount,boimupcount, by='Family') %>% full_join(.,cmimcount,by='Family') %>% full_join(.,cmimupcount,by='Family') %>% full_join(.,ncimcount,by='Family') %>% full_join(.,rimupcount,by='Family') %>% full_join(.,simupcount,by='Family') %>% full_join(.,wimupcount,by='Family')
Refinetable<- subset(Sumtable, !(boup==0 & cmup==0 &rup==0 & sup==0 & wup==0))
write.csv(Refinetable,"RegulateSum.csv",sep="\t",quote=F)

# Manually check

allgenecount<- read.csv("RegulateSum.csv",sep=",",header=T)
pdf("Regulated_immune_species.pdf",35,15)
par(mfrow = c(5, 1))
barplot(allgenecount$wup,col=c(rep("black",5), rep("green3",5),rep("red",8)),cex.names=1.5,las=1,tck=0.1,ylim=c(0,10),ylab="Worker",cex=1.5)
barplot(allgenecount$sup,col=c(rep("black",5), rep("green3",5),rep("red",8)),cex.names=1.5,las=1,tck=0.1,ylim=c(0,10),ylab="Soldier",cex=1.5)
barplot(allgenecount$rup,col=c(rep("black",5), rep("green3",5),rep("red",8)),cex.names=1.5,las=1,tck=0.1,ylim=c(0,10),ylab="Reproductives",cex=1.5)
barplot(allgenecount$cmup,col=c(rep("black",5), rep("green3",5),rep("red",8)),cex.names=1.5,las=1,tck=0.1,ylim=c(0,10),ylab="Cm",cex=1.5)
barplot(allgenecount$boup,names.arg=allgenecount$Family,col=c(rep("black",5), rep("green3",5),rep("red",8)),cex.names=1.5,las=1,tck=0.1,ylim=c(0,10),ylab="Bo",cex=1.5)
dev.off()