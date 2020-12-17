library(dplyr)
library(tidyr)



allgene<-read.csv("all_db6_immune_checkid_family.csv",sep="\t",head=F,stringsAsFactors=F) %>% subset(V2=="0"|V2 =="1") %>% tidyr::separate(V7,c("DN","Cluster"),"_c",remove=F)
allgene<- allgene[grep("quant", allgene$V1,invert=T),]

if(FALSE){
for (spec in unique(allgene$V1)){
genetable<-dplyr::filter(allgene,V1==spec)

#genetable<-read.csv("test0",sep="\t",head=F) %>% tidyr::separate(V7,c("DN","Cluster"),"_c",remove=F)'

DNreplicate<- dplyr::filter(genetable, DN%in%genetable[duplicated(genetable[c(8)]),8]) 

nrgene<- read.csv(paste("/media/shulinhe/DATA/QuanL/NR_diamond/",spec,".match",sep=""),sep="\t",head=F) %>% dplyr::filter(V1%in%DNreplicate$V7) %>% dplyr::filter(V13<=0.00001)

#top_n(n=5,evalue), slice(1:5),do(head(., n=5))"""




final<-data.frame("V1"="","no_overlap_gene"="","V2"="","V9"=0,"V10"=0,"Vs"=0,"Ve"=0)

for (i in unique(DNreplicate[,c("DN")])){

	reiso<- dplyr::filter(DNreplicate,DN==i)

	rdn<- dplyr::filter(nrgene,V1%in%reiso$V7) %>% group_by(V1) %>% dplyr::filter(row_number()<=10L) %>% arrange(V2,V9) %>% as.data.frame

	rdntarget<- dplyr::filter(rdn,V2%in%rdn[duplicated(rdn$V2),2])

	for (j in unique(rdntarget$V2)){
		rdnstarget<- dplyr::filter(rdn, V2==j)
		h<- 1
		while (h<nrow(rdnstarget)){
			g=1
			while(g<=(nrow(rdnstarget)-h)){
			if (rdnstarget[h,c("V10")]<(rdnstarget[h+g,c("V9")]+8)){
				special<-data.frame(V1=as.character(rdnstarget[h,1]),no_overlap_gene=as.character(rdnstarget[h+g,1]),V2=as.character(j),rdnstarget[h,9:10],Vs=rdnstarget[h+g,9],Ve=rdnstarget[h+g,10])
				final<- rbind(final,special)
				break
			} else {
				g=g+1
				}}
			h=h+1
		}

	}
	final<-rbind(final,data.frame("V1"="","no_overlap_gene"="","V2"="","V9"="","V10"="","Vs"="","Ve"=""))	
}
final$Species<- spec
print (length(unique(final$V1)))
write.table(final, "all_db6_immune_checkid_family_nr_refine.csv",quote=F,sep="\t",append=TRUE,na="")
}
}


for (spec in unique(allgene$V1)){
genetable<-dplyr::filter(allgene,V1==spec)
#DNreplicate<- dplyr::filter(genetable, DN%in%genetable[duplicated(genetable[c(8)]),8]) 

nrgene<- read.csv(paste("/media/shulinhe/DATA/QuanL/NR_diamond/",spec,".match",sep=""),sep="\t",head=F) %>% dplyr::filter(V1%in%genetable$V7)
# %>% dplyr::filter(V13<=0.00001)

#top_n(n=5,evalue), slice(1:5),do(head(., n=5))"""

final<-data.frame("V1"="","no_overlap_gene"="","V2"="","V9"=0,"V10"=0,"Vs"=0,"Ve"=0)

for (i in unique(genetable[,"V4"])){

	reiso<- dplyr::filter(genetable,V4==i)

	rdn<- dplyr::filter(nrgene,V1%in%reiso$V7) %>% group_by(V1) %>% dplyr::filter(row_number()<=6L) %>% arrange(V2,V9) %>% as.data.frame

	rdntarget<- dplyr::filter(rdn,V2%in%rdn[duplicated(rdn$V2),2])

	for (j in unique(rdntarget$V2)){
		rdnstarget<- dplyr::filter(rdn, V2==j)
		h<- 1
		while (h<nrow(rdnstarget)){
			g=1
			while(g<=(nrow(rdnstarget)-h)){
			if (rdnstarget[h,c("V10")]<(rdnstarget[h+g,c("V9")]+8)){
				special<-data.frame(V1=as.character(rdnstarget[h,1]),no_overlap_gene=as.character(rdnstarget[h+g,1]),V2=as.character(j),rdnstarget[h,9:10],Vs=rdnstarget[h+g,9],Ve=rdnstarget[h+g,10])
				final<- rbind(final,special)
				break
			} else {
				g=g+1
				}}
			h=h+1
		}

	}
	final<-rbind(final,data.frame("V1"="","no_overlap_gene"="","V2"="","V9"="","V10"="","Vs"="","Ve"=""))	
}
final$Species<- spec
print (length(unique(final$V1)))
write.table(final, "linshi",quote=F,sep="\t",append=TRUE,na="")
}

