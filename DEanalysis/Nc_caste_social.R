#!/usr/bin/Rscript

library(tximport)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(tibble)
library(DESeq2)
library(apeglm)
library(pheatmap)

dir <- "/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/"
tx2gene <- read_tsv("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/tx2gene.tsv")

library <- c("CR1_S5","CR2_S11","CR3_S17","CS1_S3","CS2_S9","CS3_S15","CW1_S1","CW2_S7","CW3_S13","TR1_S6","TR2_S12","TR3_S18","TS1_S4","TS2_S10","TS3_S16","TW1_S2","TW2_S8","TW3_S14")
files <- file.path(dir, library, "quant.sf")
names(files) <- library
all(file.exists(files))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("treatment", "rest"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^C", "control_", sampleTable$treatment)
sampleTable$treatment <- gsub("^T", "treatment_", sampleTable$treatment)
sampleTable <- separate(sampleTable, treatment, into = c("treatment", "caste"), sep = "_", extra = "drop")
sampleTable$caste <- gsub("^R.", "reproductive", sampleTable$caste)
sampleTable$caste <- gsub("^S.", "soldier", sampleTable$caste)
sampleTable$caste <- gsub("^W.", "worker", sampleTable$caste)
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, -library, -rest)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment+caste)
bpgene<- read.csv("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_full_tax_nr_metazoa",header=F,sep="\t")

dds<- dds[rownames(dds) %in% bpgene$V1,]
dds <- DESeq(dds, parallel = TRUE, fitType = "local")
#dds_lrt<-DESeq(dds, test="LRT",reduced=~treatment, fitType="local",parallel=TRUE)

rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
pdf("pca.pdf")
plot(pca)
dev.off()
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
pcax$caste <- sampleTable$caste
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>% .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment,shape=caste)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment,shape=caste)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment,shape=caste)) + geom_point(size=5) + xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("social_caste_PCA.pdf", 21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

Immune<-read.table("/media/shulinhe/DATA/Full_analysis/immune/nc_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(Immune)<- c("abc","Gene_ID","Family","hmm1","hmm2","iso","type","len","val","protid","blastname","blastiden","fullname","tax","pfamd","pfamdname","pfame")
pca <- prcomp(t(rldf[rownames(rldf)%in%Immune$Gene_ID,]))
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
pcax$caste <- sampleTable$caste
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>% .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment,shape=caste)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment,shape=caste)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment,shape=caste)) + geom_point(size=5) + xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("social_caste_immune_PCA.pdf", 21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()
allimmuneheat<- assay(rld)[as.character(Immune$Gene_ID),c(1:3,10:12,4:6,13:15,7:9,16:18)]
allimmuneid<- order(rowMeans(allimmuneheat[,c(1:6)]),decreasing=T)
pdf("social_caste_immune_heatmap.pdf")
pheatmap(allimmuneheat[allimmuneid,],annotation_col=sampleTable,cluster_row=F,show_rownames = F)
dev.off()

#################caste comparasion#########################
SR_MA <- lfcShrink(dds, coef=3, type="apeglm")
WR_MA <- lfcShrink(dds, coef=4, type="apeglm")
dds$caste<- relevel(dds$caste, "worker")
dds <- nbinomWaldTest(dds)
WS_MA <- lfcShrink(dds, coef=, type="normal")
tiff("caste_MAplot.tiff",units="in",8,4,re=600)
par(mfrow=c(1, 3))
plotMA(WR_MA, ylim=c(-8,8),main="Workers vs Reproductives",cex=0.55)
plotMA(SR_MA, ylim=c(-8,8),main="Soldiers vs Reproductives",cex=0.55)
plotMA(SW_MA, ylim=c(-8,8),main="Soldiers vs Workers",cex=0.55)
dev.off()

RS_res_sig<- results(dds, alpha = 0.01,contrast=c("caste","reproductive","soldier"),tidy=TRUE)
RvsS_sig_up <- dplyr::filter(as.data.frame(RS_res_sig), padj < 0.01 & log2FoldChange > 2)
write_csv(RvsS_sig_up,"./RvsS_sig_up.csv")
RvsS_sig_down <- dplyr::filter(as.data.frame(RS_res_sig), padj < 0.01 & log2FoldChange < -2)
write_csv(RvsS_sig_down,"./RvsS_sig_down.csv")

RW_res_sig<- results(dds, alpha = 0.01,contrast=c("caste","reproductive","worker"),tidy=TRUE) #%>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
RvsW_sig_up <- dplyr::filter(as.data.frame(RW_res_sig), padj < 0.01 & log2FoldChange > 2)
write_csv(RvsW_sig_up,"./RvsW_sig_up.csv")
RvsW_sig_down <- dplyr::filter(as.data.frame(RW_res_sig), padj < 0.01 & log2FoldChange < -2)
write_csv(RvsW_sig_down,"./RvsW_sig_down.csv")

WS_res_sig<- results(dds, alpha = 0.01,contrast=c("caste","worker","soldier"),tidy=TRUE) #%>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
WvsS_sig_up <- dplyr::filter(as.data.frame(WS_res_sig), padj < 0.01 & log2FoldChange > 2)
write_csv(WvsS_sig_up,"./WvsS_sig_up.csv")
WvsS_sig_down <- dplyr::filter(as.data.frame(WS_res_sig), padj < 0.01 & log2FoldChange < -2)
write_csv(WvsS_sig_down,"./WvsS_sig_down.csv")

RvsWS_sig_up<- RvsS_sig_up[RvsS_sig_up$row %in%RvsW_sig_up$row,]
write_csv(RvsWS_sig_up,"./RvsWS_sig_up.csv")
RvsWS_sig_down<- RvsS_sig_down[RvsS_sig_down$row %in%RvsW_sig_down$row,]
write_csv(RvsWS_sig_down,"./RvsWS_sig_down.csv")

SvsRW_sig_up<- RvsS_sig_down[RvsS_sig_down$row %in%WvsS_sig_down$row,]
write_csv(SvsRW_sig_up,"./SvsRW_sig_up.csv")
SvsRW_sig_down<- RvsS_sig_up[RvsS_sig_up$row %in%WvsS_sig_up$row,]
write_csv(SvsRW_sig_down,"./SvsRW_sig_down.csv")

WvsRS_sig_up<- RvsW_sig_down[RvsW_sig_down$row %in%WvsS_sig_up$row,]
write_csv(WvsRS_sig_up,"./WvsRS_sig_up.csv")
WvsRS_sig_down<- RvsW_sig_up[RvsW_sig_up$row %in%WvsS_sig_down$row,]
write_csv(WvsRS_sig_down,"./WvsRS_sig_down.csv")

RWupim<- RvsW_sig_up[RvsW_sig_up$row %in% Immune$Gene_ID,]
RSupim<- RvsS_sig_up[RvsS_sig_up$row %in% Immune$Gene_ID,]
RWdownim<- RvsW_sig_down[RvsW_sig_down$row %in% Immune$Gene_ID,]
RSdownim<- RvsS_sig_down[RvsS_sig_down$row %in% Immune$Gene_ID,]
WSupim<- WvsS_sig_up[WvsS_sig_up$row %in% Immune$Gene_ID, ]
WSdownim<- WvsS_sig_down[WvsS_sig_down$row %in% Immune$Gene_ID, ]

refe<- assay(rld)["TRINITY_DN139618_c0_g1",]
allcastesImmune_regulate<- c(as.character(RWupim[,c("row")]),as.character(RSupim[,c("row")]),as.character(RWdownim[,c("row")]),as.character(RSdownim[,c("row")]),as.character(WSupim[,c("row")]),as.character(WSdownim[,c("row")]))
allImmune_caste_regulate<- unique(assay(rld)[allcastesImmune_regulate,])
casteorder<- order(rowSums(allImmune_caste_regulate[,c(1:3,10:12)]), decreasing=T)
totalcaste<- rbind(refe,allImmune_caste_regulate[casteorder,])
totalcaste<- totalcaste[,c(1:3,10:12,4:6,13:15,7:9,16:18)]
library(pheatmap)
pdf("social_caste_comparasion_immune_heatmap.pdf",10,14)
pheatmap(totalcaste,cluster_rows=F, annotation_col=sampleTable,cluster_cols=F)
dev.off()

if (FALSE){ 
allcastesImmune_regulate<- c(as.character(RSupim[,c("row")]),as.character(RSdownim[,c("row")]),as.character(WSupim[,c("row")]),as.character(WSdownim[,c("row")]))
allImmune_caste_regulate<- assay(rld)[allcastesImmune_regulate,]
totalcaste<- rbind(refe,allImmune_caste_regulate)
totalcaste<- totalcaste[,c(1:3,10:12,4:6,13:15,7:9,16:18)]
library(pheatmap)
pdf("social_caste_S_immune_heatmap.pdf",10,14)
pheatmap(totalcaste,cluster_rows=F, annotation_col=sampleTable,cluster_cols=F)
dev.off()

allcastesImmune_regulate<- c(as.character(RWupim[,c("row")]),as.character(RWdownim[,c("row")]),as.character(WSupim[,c("row")]),as.character(WSdownim[,c("row")]))
allImmune_caste_regulate<- assay(rld)[allcastesImmune_regulate,]
totalcaste<- rbind(refe,allImmune_caste_regulate)
totalcaste<- totalcaste[,c(1:3,10:12,4:6,13:15,7:9,16:18)]
library(pheatmap)
pdf("social_caste_W_immune_heatmap.pdf",10,14)
pheatmap(totalcaste,cluster_rows=F, annotation_col=sampleTable,cluster_cols=F)
dev.off()
}


RWupim$note<- "RupW"
RSupim$note<- "RupS"
RWdownim$note<- "RdownW"
RSdownim$note<- "RdownS"
WSupim$note<- "WupS"
WSdownim$note<- "WdownS"

write_csv(rbind(RWupim,RSupim,RWdownim,RSdownim,WSupim,WSdownim), "./social_caste_comparasion_immune.csv")

################treatment comparasion######################
TC_MA <- lfcShrink(dds, coef=2, type="apeglm")
png("TC_MAplot")
plotMA(TC_MA, ylim=c(-7,7),main="treatment_control_total_nestmate",cex=0.8)
dev.off()
TC_res_sig<- results(dds, alpha = 0.05,contrast=c("treatment","treatment","control"),tidy=TRUE) #%>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
TC_sig_up <- dplyr::filter(as.data.frame(TC_res_sig), padj < 0.05 & log2FoldChange > 1)
write_csv(TC_sig_up,"./TC_sig_up.csv")
TC_sig_down <- dplyr::filter(as.data.frame(TC_res_sig), padj < 0.05 & log2FoldChange < -1)
write_csv(TC_sig_down,"./TC_sig_down.csv")

###############treatment*caste comparasion#################

sampleTable_i <- unite(sampleTable, treatment_caste, c(treatment, caste), remove = FALSE)
ddi <- DESeqDataSetFromTximport(txi, sampleTable_i, ~ treatment_caste)
ddi<- ddi[rownames(ddi) %in% bpgene$V1,]
ddi <- DESeq(ddi, parallel = TRUE, fitType = "local")

#TC_R_MA <- lfcShrink(ddi, contrast=c("treatment_caste","treatment_reproductive","control_reproductive"), type="normal")
TC_R_MA <- lfcShrink(ddi, coef=4, lfcThreshold=1, type="apeglm")

ddi$treatment_caste<- relevel(ddi$treatment_caste, "control_soldier")
ddi <- nbinomWaldTest(ddi)
#TC_S_MA <- lfcShrink(ddi, contrast=c("treatment_caste","treatment_soldier","control_soldier"), type="normal")
TC_S_MA <- lfcShrink(ddi, coef=5, type="apeglm")

ddi$treatment_caste<- relevel(ddi$treatment_caste, "control_worker")
ddi <- nbinomWaldTest(ddi)
#TC_W_MA <- lfcShrink(ddi, contrast=c("treatment_caste","treatment_worker","control_worker"), type="normal")
TC_W_MA <- lfcShrink(ddi, coef=6, type="apeglm")
tiff("caste_R_treatment_MAplot.tiff",units="in",width=4,height=4,re=600)
plotMA(TC_R_MA, ylim=c(-6,6),main="Reproductive",cex=0.55)
dev.off()
tiff("caste_S_treatment_MAplot.tiff",units="in",width=4,height=4,re=600)
plotMA(TC_S_MA, ylim=c(-6,6),main="Soldier",cex=0.55)
dev.off()
tiff("caste_W_treatment_MAplot.tiff",units="in",width=4,height=4,re=600)
plotMA(TC_W_MA, ylim=c(-7,7),main="Worker",cex=0.55)
dev.off()

TC_R_res_sig<- results(ddi, alpha = 0.05,contrast=c("treatment_caste","treatment_reproductive","control_reproductive"),tidy=TRUE) # %>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
TC_R_sig_up <- dplyr::filter(as.data.frame(TC_R_res_sig), padj < 0.05 & log2FoldChange > 1)
write_csv(TC_R_sig_up,"./TC_R_sig_up.csv")
TC_R_sig_down <- dplyr::filter(as.data.frame(TC_R_res_sig), padj < 0.05 & log2FoldChange < -1)
write_csv(TC_R_sig_down,"./TC_R_sig_down.csv")

TC_S_res_sig<- results(ddi, alpha = 0.05,contrast=c("treatment_caste","treatment_soldier","control_soldier"),tidy=TRUE) # %>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
TC_S_sig_up <- dplyr::filter(as.data.frame(TC_S_res_sig), padj < 0.05 & log2FoldChange > 1)
write_csv(TC_S_sig_up,"./TC_S_sig_up.csv")
TC_S_sig_down <- dplyr::filter(as.data.frame(TC_S_res_sig), padj < 0.05 & log2FoldChange < -1)
write_csv(TC_S_sig_down,"./TC_S_sig_down.csv")

TC_W_res_sig<- results(ddi, alpha = 0.05,contrast=c("treatment_caste","treatment_worker","control_worker"),tidy=TRUE) #%>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
TC_W_sig_up <- dplyr::filter(as.data.frame(TC_W_res_sig), padj < 0.05 & log2FoldChange > 1)
write_csv(TC_W_sig_up,"./TC_W_sig_up.csv")
TC_W_sig_down <- dplyr::filter(as.data.frame(TC_W_res_sig), padj < 0.05 & log2FoldChange < -1)
write_csv(TC_W_sig_down,"./TC_W_sig_down.csv")

####immune gene expression plot####

rImmune_up<- TC_R_sig_up[TC_R_sig_up$row %in% Immune$Gene_ID,]
rImmune_down<- TC_R_sig_down[TC_R_sig_down$row %in% Immune$Gene_ID,]
sImmune_up<- TC_S_sig_up[TC_S_sig_up$row %in% Immune$Gene_ID,]
sImmune_down<- TC_S_sig_down[TC_S_sig_down$row %in% Immune$Gene_ID,]
wImmune_up<- TC_W_sig_up[TC_W_sig_up$row %in% Immune$Gene_ID,]
wImmune_down<- TC_W_sig_down[TC_W_sig_down$row %in% Immune$Gene_ID,]
allsImmune_regulate<- c(as.character(rImmune_up[,c("row")]),as.character(sImmune_up[,c("row")]),as.character(wImmune_up[,c("row")]),as.character(rImmune_down[,c("row")]),as.character(sImmune_down[,c("row")]),as.character(wImmune_down[,c("row")]))
allImmune_regulate<- assay(rld)[allsImmune_regulate,]
refe<- assay(rld)["TRINITY_DN139618_c0_g1",]
total<- rbind(refe,allImmune_regulate)
total<- total[,c(1:3,10:12,4:6,13:15,7:9,16:18)]
library(pheatmap)
pdf("social_nestmate_regulate_immune_heatmap.pdf",10,14)
pheatmap(total,cluster_rows=F, annotation_col=sampleTable,cluster_cols=F)
dev.off()



