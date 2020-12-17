#!/usr/bin/Rscript

library(tximport)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(tibble)
library(DESeq2)
library(apeglm)
library(IHW)

dir <- "/media/shulinhe/DATA/Full_analysis/Bo_full/bo_salmon/"
#############################nestmates################################

#if (FALSE){ 
run <- c("Cn1_S7", "Cn2_S8", "Cn3_S9", "Tn1_S10","Tn2_S11","Tn3_S12")
files <- file.path(dir, run, "quant.sf")
names(files) <- run
all(file.exists(files))
tx2gene <- read_tsv("/media/shulinhe/DATA/Full_analysis/Bo_full/bo_tx2gene.tsv")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
rownames(sampleTable) <- sampleTable$library
sampleTable <- separate(sampleTable, library, into = c("treatment", "libID"), sep = "_", extra = "drop")
sampleTable$treatment <- gsub("^C..$", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^T..$", "treatment", sampleTable$treatment)
sampleTable <- dplyr::select(sampleTable, -libID)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment)

bpgene<- read.csv("/media/shulinhe/DATA/Full_analysis/Bo_full/bo_full_tax_nr_metazoa",header=F,sep="\t")
dds<- dds[rownames(dds) %in% bpgene$V1,]

dds <- DESeq(dds, parallel = TRUE, fitType = "local")
resLFC <- lfcShrink(dds, coef=2,type="apeglm")
tiff("bo_n_MAplot.tiff",units="in",res=600,width=4,height=4)
plotMA(resLFC, ylim=c(-4,4),cex=0.55)
dev.off()

rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
pdf("bo_n_pca.pdf")
plot(pca)
dev.off()

pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>%
  .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("bo_n_PCA.pdf",21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

#res <- results(dds, alpha = 0.01,filterFun="ihw") %>% subset(padj < 0.01 & abs(log2FoldChange) > 1)
#res <- results(dds, alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")#
res_tvsc <- results(dds, alpha = 0.05, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
res_TvsC_sig_up <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.05 & log2FoldChange > 1)
write_csv(res_TvsC_sig_up,"./bon_res_TvsC_sig_up.csv")
res_TvsC_sig_down <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.05 & log2FoldChange < -1)
write_csv(res_TvsC_sig_down,"./bon_res_TvsC_sig_down.csv")
#}
##immune gene plot
#if (TRUE){ 
Immune<-read.table("/media/shulinhe/DATA/Full_analysis/immune/bo_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(Immune)<- c("clust","Gene_ID","Family","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
Immune_up<- res_TvsC_sig_up[res_TvsC_sig_up$row %in% Immune$Gene_ID,]
Immune_down<- res_TvsC_sig_down[res_TvsC_sig_down$row %in% Immune$Gene_ID,]
library(pheatmap)

immupe<- assay(rld)[as.character(Immune_up[,c("row")]),]
immupeod<- order(rowSums(immupe[,c(4,5,6)]),decreasing=T)
if (nrow(Immune_down)>1){
immude<- assay(rld)[as.character(Immune_down[,c("row")]),]
immude<- immude[order(rowSums(immude[,c(4,5,6)]),decreasing=T),]
} else {
immude<- assay(rld)[as.character(Immune_down[,c("row")]),] %>% data.frame %>% t
rownames(immude)<- as.character(Immune_down[,c("row")])
}


refe<- assay(rld)["TRINITY_DN59362_c0_g4",]

total<- rbind(refe,immupe[immupeod,], immude)
#totalname<- Immune[Immune$Gene_ID %in% rownames(total),]
#totalname<- dplyr::select(totalname, "Gene_ID","Blastp_target_full_name")
#total1<- cbind(total,Immune[match(rownames(total),Immune$Gene_ID),][,c(2,8)])
pdf("bon_regulate_immune_heatmap.pdf",10,14)
pheatmap(total,cluster_rows=F, annotation_col=sampleTable)
dev.off()
#}
################################injected individuals in social############################################
run <- c("Ci1_S1", "Ci2_S2", "Ci3_S3", "Ti1_S4","Ti2_S5","Ti3_S6")
files <- file.path(dir, run, "quant.sf")
names(files) <- run
all(file.exists(files))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
rownames(sampleTable) <- sampleTable$library
sampleTable <- separate(sampleTable, library, into = c("treatment", "libID"), sep = "_", extra = "drop")
sampleTable$treatment <- gsub("^C..$", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^T..$", "treatment", sampleTable$treatment)
sampleTable <- dplyr::select(sampleTable, -libID)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment)

bpgene<- read.csv("/media/shulinhe/DATA/Full_analysis/Bo_full/bo_full_tax_nr_metazoa",header=F,sep="\t")
dds<- dds[rownames(dds) %in% bpgene$V1,]

dds <- DESeq(dds, parallel = TRUE, fitType = "local")
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
png("bo_i_MAplot.png")
plotMA(resLFC, ylim=c(-7,7))
dev.off()

rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
pdf("bo_i_pca.pdf")
plot(pca)
dev.off()

pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>%
  .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("bo_i_PCA.pdf",21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

#res <- results(dds, alpha = 0.01) %>% subset(padj < 0.01 & abs(log2FoldChange) > 1)

#res <- results(dds, alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")#

res_tvsc <- results(dds, alpha = 0.01, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
res_TvsC_sig_up <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.01 & log2FoldChange > 1)
write_csv(res_TvsC_sig_up,"./boi_res_TvsC_sig_up.csv")
res_TvsC_sig_down <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.01 & log2FoldChange < -1)
write_csv(res_TvsC_sig_down,"./boi_res_TvsC_sig_down.csv")

#if (TRUE){
Immune<-read.table("../../immune/bo_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(Immune)<- c("clust","Gene_ID","Family","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
Immune_up<- res_TvsC_sig_up[res_TvsC_sig_up$row %in% Immune$Gene_ID,]
Immune_down<- res_TvsC_sig_down[res_TvsC_sig_down$row %in% Immune$Gene_ID,]
library(pheatmap)
library(pheatmap)

immupe<- assay(rld)[as.character(Immune_up[,c("row")]),]
immupeod<- order(rowSums(immupe[,c(4,5,6)]),decreasing=T)
if (nrow(Immune_down)>1){
immude<- assay(rld)[as.character(Immune_down[,c("row")]),]
immude<- immude[order(rowSums(immude[,c(4,5,6)]),decreasing=T),]
} else {
immude<- assay(rld)[as.character(Immune_down[,c("row")]),] %>% data.frame %>% t
rownames(immude)<- as.character(Immune_down[,c("row")])
}


refe<- assay(rld)["TRINITY_DN59362_c0_g4",]

total<- rbind(refe,immupe[immupeod,], immude)
#totalname<- Immune[Immune$Gene_ID %in% rownames(total),]
#totalname<- dplyr::select(totalname, "Gene_ID","Blastp_target_full_name")
#total1<- cbind(total,Immune[match(rownames(total),Immune$Gene_ID),][,c(2,8)])
pdf("boi_regulate_immune_heatmap.pdf",10,14)
pheatmap(total,cluster_rows=F, annotation_col=sampleTable)
dev.off()
#}



############################injected individuals##################################

run <- c("BoC1_S1", "BoC2_S2", "BoT1_S3", "BoT2_S4")
files <- file.path(dir, run, "quant.sf")
names(files) <- run
all(file.exists(files))
tx2gene <- read_tsv("bo_tx2gene.tsv")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
#'''this derives the sample names from the count dataframe - safer than doing it manually'''
sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("treatment", "species"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^BoC.$", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^BoT.$", "treatment", sampleTable$treatment)
rownames(sampleTable) <- sampleTable$library
#'''drop variables that are not be fitted'''
sampleTable <- dplyr::select(sampleTable, -library, -species)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment)

#'''filter bacteria and protist genes'''
dds<- dds[rownames(dds) %in% bpgene$V1,]

dds <- DESeq(dds, parallel = TRUE, fitType = "local")
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
png("bo_ii_MAplot.png")
plotMA(resLFC, ylim=c(-7,7))
dev.off()

rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
#'''scree plot of PCs'''
pdf("bo_ii_pca.pdf")
plot(pca)
dev.off()
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>%
  .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("bo_ii_pca123.pdf",21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

#res <- results(dds, alpha = 0.05) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
#res <- results(dds, alpha = 0.05, lfcThreshold=1, altHypothesis="greaterAbs")#
res_tvsc <- results(dds, alpha = 0.01, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
res_TvsC_sig_up <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.01 & log2FoldChange > 1)
write_csv(res_TvsC_sig_up,"./bo_ii_TvsC_sig_up.csv")
res_TvsC_sig_down <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.01 & log2FoldChange < -1)
write_csv(res_TvsC_sig_down,"./bo_ii_TvsC_sig_down.csv")

#if (TRUE){
Immune<-read.table("../../immune/bo_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(Immune)<- c("clust","Gene_ID","Family","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
Immune_up<- res_TvsC_sig_up[res_TvsC_sig_up$row %in% Immune$Gene_ID,]
Immune_down<- res_TvsC_sig_down[res_TvsC_sig_down$row %in% Immune$Gene_ID,]
library(pheatmap)
library(pheatmap)

immupe<- assay(rld)[as.character(Immune_up[,c("row")]),]
immupeod<- order(rowSums(immupe[,c(3,4)]),decreasing=T)
if (nrow(Immune_down)>1){
immude<- assay(rld)[as.character(Immune_down[,c("row")]),]
immude<- immude[order(rowSums(immude[,c(3,4)]),decreasing=T),]
} else {
immude<- assay(rld)[as.character(Immune_down[,c("row")]),] %>% data.frame %>% t
rownames(immude)<- as.character(Immune_down[,c("row")])
}


refe<- assay(rld)["TRINITY_DN59362_c0_g4",]

total<- rbind(refe,immupe[immupeod,], immude)
#totalname<- Immune[Immune$Gene_ID %in% rownames(total),]
#totalname<- dplyr::select(totalname, "Gene_ID","Blastp_target_full_name")
#total1<- cbind(total,Immune[match(rownames(total),Immune$Gene_ID),][,c(2,8)])
pdf("bo_ii_regulate_immune_heatmap.pdf",10,14)
pheatmap(total,cluster_rows=F, annotation_col=sampleTable)
dev.off()
#}
