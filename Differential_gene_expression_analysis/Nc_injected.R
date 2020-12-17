#!/usr/bin/Rscript

library(tximport)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(tibble)
library(DESeq2)
library(apeglm)

dir <- "/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/"
tx2gene <- read_tsv("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/tx2gene.tsv")

#######################Injected workers in social condition################################
library <- c("C1_S1","C2_S2","C3_S3","T1_S4","T2_S5","T3_S6")
files <- file.path(dir, library, "quant.sf")
names(files) <- library
all(file.exists(files))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("treatment", "rest"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^C.", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^T.", "treatment", sampleTable$treatment)
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, -library, -rest)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment)
bpgene<- read.csv("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_full_tax_nr_metazoa",header=F,sep="\t")
dds<- dds[rownames(dds) %in% bpgene$V1,]
dds <- DESeq(dds, parallel = TRUE, fitType = "local")

rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
pdf("Injected_social_pca.pdf")
plot(pca)
dev.off()
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>% .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("Injected_social_PCA.pdf", 21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

I_S_MA<- lfcShrink(dds, coef=2, type="apeglm")
png("I_S_MAplot")
plotMA(I_S_MA, ylim=c(-7,7),cex=0.8,main="Injected workers in social")
dev.off()

#I_S_res_sig<- results(dds, alpha = 0.01, tidy=TRUE) %>% subset(padj < 0.01 & abs(log2FoldChange) > 1)
IS_tvsc <- results(dds, alpha = 0.01, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
IS_TvsC_sig_up <- dplyr::filter(as.data.frame(IS_tvsc), padj < 0.01 & log2FoldChange > 1)
write_csv(IS_TvsC_sig_up,"./IS_TvsC_sig_up.csv")
IS_TvsC_sig_down <- dplyr::filter(as.data.frame(IS_tvsc), padj < 0.01 & log2FoldChange < -1)
write_csv(IS_TvsC_sig_down,"./IS_TvsC_sig_down.csv")

Immune<-read.table("/media/shulinhe/DATA/Full_analysis/immune/nc_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(Immune)<- c("clust","Gene_ID","Family","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
wsImmune_up<- IS_TvsC_sig_up[IS_TvsC_sig_up$row %in% Immune$Gene_ID,]
wsImmune_down<- IS_TvsC_sig_down[IS_TvsC_sig_down$row %in% Immune$Gene_ID,]
wsImmune_regulate<- c(as.character(wsImmune_up[,c("row")]),as.character(wsImmune_up[,c("row")]))
allsImmune_regulate<- assay(rld)[wsImmune_regulate,]
refe<- assay(rld)["TRINITY_DN139618_c0_g1",]
totals<- rbind(refe,allsImmune_regulate)
library(pheatmap)
pdf("nc_worker_si_regulate_immune_heatmap.pdf",10,14)
pheatmap(totals,cluster_rows=F, annotation_col=sampleTable)
dev.off()

#######################Injected workers in individual condition################################
#if (FALSE){
wrun <- c("nc_cp_1_S3","nc_cp_2_S11","nc_tp_1_S4","nc_tp_2_S12")
wfiles <- file.path(dir, wrun, "quant.sf")
names(wfiles) <- wrun
all(file.exists(wfiles))

wtxi <- tximport(wfiles, type = "salmon", tx2gene = tx2gene)
sampleTable <- as.data.frame(colnames(wtxi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("species", "treatment"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^cp$", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^tp$", "treatment", sampleTable$treatment)
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, -library, -species)
wdds <- DESeqDataSetFromTximport(wtxi, sampleTable, ~ treatment)


wdds<- wdds[rownames(wdds) %in% bpgene$V1,]
wdds <- DESeq(wdds, parallel = TRUE, fitType = "local")

wresLFC <- lfcShrink(wdds, coef=2, type="apeglm")
png("wMAplot")
plotMA(wresLFC, ylim=c(-7,7),cex=0.8,main="individual injected workers")
dev.off()

wrld <- rlog(wdds, fitType = "local")
wldf <- assay(wrld) %>% as.data.frame
pca <- prcomp(t(wldf))
pdf("wpca.pdf")
plot(pca)
dev.off()
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>% .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("wPCA.pdf", 21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

#wres <- results(wdds, alpha = 0.05) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
wres_tvsc <- results(wdds, alpha = 0.01, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
wres_TvsC_sig_up <- dplyr::filter(as.data.frame(wres_tvsc), padj < 0.01 & log2FoldChange > 1)
write_csv(wres_TvsC_sig_up,"./wres_TvsC_sig_up.csv")
wres_TvsC_sig_down <- dplyr::filter(as.data.frame(wres_tvsc), padj < 0.01 & log2FoldChange < -1)
write_csv(wres_TvsC_sig_down,"./wres_TvsC_sig_down.csv")

#######################Injected reproductives in individual condition################################

rrun <- c("nc_cr_1_S7","nc_cr_2_S15","nc_tr_1_S8","nc_tr_2_S16")
rfiles <- file.path(dir, rrun, "quant.sf")
names(rfiles) <- rrun
all(file.exists(rfiles))

rtxi <- tximport(rfiles, type = "salmon", tx2gene = tx2gene)
sampleTable <- as.data.frame(colnames(rtxi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("species", "treatment"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^cr$", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^tr$", "treatment", sampleTable$treatment)
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, -library, -species)
rdds <- DESeqDataSetFromTximport(rtxi, sampleTable, ~ treatment)
rdds<- rdds[rownames(rdds) %in% bpgene$V1,]
rdds <- DESeq(rdds, parallel = TRUE, fitType = "local")

rresLFC <- lfcShrink(rdds, coef=2,type="apeglm")
png("rMAplot")
plotMA(rresLFC, ylim=c(-7,7),cex=0.8,main="individual injected reproductives")
dev.off()

rrld <- rlog(rdds, fitType = "local")
rldf <- assay(rrld) %>% as.data.frame
pca <- prcomp(t(rldf))
pdf("rpca.pdf")
plot(pca)
dev.off()
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>% .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("rPCA.pdf", 21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

#rres <- results(rdds, alpha = 0.05) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
rres_tvsc <- results(rdds, alpha = 0.01, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
rres_TvsC_sig_up <- dplyr::filter(as.data.frame(rres_tvsc), padj < 0.01 & log2FoldChange > 1)
write_csv(rres_TvsC_sig_up,"./rres_TvsC_sig_up.csv")
rres_TvsC_sig_down <- dplyr::filter(as.data.frame(rres_tvsc), padj < 0.01 & log2FoldChange < -1)
write_csv(rres_TvsC_sig_down,"./rres_TvsC_sig_down.csv")


#######################Injected soldiers in individual condition################################
srun <- c("nc_cs_1_S5","nc_cs_2_S13","nc_ts_1_S6","nc_ts_2_S14")
sfiles <- file.path(dir, srun, "quant.sf")
names(sfiles) <- srun
all(file.exists(sfiles))
stxi <- tximport(sfiles, type = "salmon", tx2gene = tx2gene)
sampleTable <- as.data.frame(colnames(stxi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("species", "treatment"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^cs$", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^ts$", "treatment", sampleTable$treatment)
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, -library, -species)
sdds <- DESeqDataSetFromTximport(stxi, sampleTable, ~ treatment)

sdds<- sdds[rownames(sdds) %in% bpgene$V1,]
sdds <- DESeq(sdds, parallel = TRUE, fitType = "local")

sresLFC <- lfcShrink(sdds, coef=2,type="apeglm")
png("sMAplot")
plotMA(sresLFC, ylim=c(-7,7),cex=0.8, main="individual injected soldiers")
dev.off()

srld <- rlog(sdds, fitType = "local")
sldf <- assay(srld) %>% as.data.frame
pca <- prcomp(t(sldf))
pdf("spca.pdf")
plot(pca)
dev.off()
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>% .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) + xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) + ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) + theme(legend.justification=c(1,0), legend.position=c(1,0))
pdf("sPCA.pdf", 21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

#sres <- results(sdds, alpha = 0.05) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1)
sres_tvsc <- results(sdds, alpha = 0.01, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
sres_TvsC_sig_up <- dplyr::filter(as.data.frame(sres_tvsc), padj < 0.01 & log2FoldChange > 1)
write_csv(sres_TvsC_sig_up,"./sres_TvsC_sig_up.csv")
sres_TvsC_sig_down <- dplyr::filter(as.data.frame(sres_tvsc), padj < 0.01 & log2FoldChange < -1)
write_csv(sres_TvsC_sig_down,"./sres_TvsC_sig_down.csv")
#}


#######################immune gene regulation################################
rm(srld,sldf, pca, rrld,rldf,wrld,wldf,sdds,rdds,wdds,rld)

run <- readLines("/media/shulinhe/DATA/Full_analysis/Nc_full/nc_salmon/nc_samples.txt")
files <- file.path(dir, run, "quant.sf")
names(files) <- run
all(file.exists(files))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("species", "treatment"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^c", "control_", sampleTable$treatment)
sampleTable$treatment <- gsub("^t", "treatment_", sampleTable$treatment)
sampleTable <- separate(sampleTable, treatment, into = c("treatment", "caste"), sep = "_")
sampleTable$caste <- gsub("^p$", "worker", sampleTable$caste)
sampleTable$caste <- gsub("^r$", "reproductive", sampleTable$caste)
sampleTable$caste <- gsub("^s$", "soldier", sampleTable$caste)
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, treatment, caste)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment * caste)
dds<- dds[rownames(dds) %in% bpgene$V1,]
dds <- DESeq(dds, parallel = TRUE, fitType = "local")
rld <- rlog(dds, fitType = "local")

rImmune_up<- rres_TvsC_sig_up[rres_TvsC_sig_up$row %in% Immune$Gene_ID,]
rImmune_down<- rres_TvsC_sig_down[rres_TvsC_sig_down$row %in% Immune$Gene_ID,]
sImmune_up<- sres_TvsC_sig_up[sres_TvsC_sig_up$row %in% Immune$Gene_ID,]
sImmune_down<- sres_TvsC_sig_down[sres_TvsC_sig_down$row %in% Immune$Gene_ID,]
wImmune_up<- wres_TvsC_sig_up[wres_TvsC_sig_up$row %in% Immune$Gene_ID,]
wImmune_down<- wres_TvsC_sig_down[wres_TvsC_sig_down$row %in% Immune$Gene_ID,]
allsImmune_regulate<- c(as.character(rImmune_up[,c("row")]),as.character(sImmune_up[,c("row")]),as.character(wImmune_up[,c("row")]),as.character(rImmune_down[,c("row")]),as.character(sImmune_down[,c("row")]),as.character(wImmune_down[,c("row")]))
allImmune_regulate<- assay(rld)[allsImmune_regulate,]
refe<- assay(rld)["TRINITY_DN139618_c0_g1",]
total<- rbind(refe,allImmune_regulate)
library(pheatmap)
pdf("nc_individual_regulate_immune_heatmap.pdf",10,14)
pheatmap(total,cluster_rows=F, annotation_col=sampleTable,cluster_cols=F)
dev.off()


if(FALSE){
immupe<- assay(rld)[as.character(Immune_up[,c("row")]),]
immude<- assay(rld)[as.character(Immune_down[,c("row")]),]
immupeod<- order(rowSums(immupe[,c(3,4)]),decreasing=T)
immudeod<- order(rowSums(immude[,c(3,4)]),decreasing=T)



total<- rbind(refe,immupe[immupeod,], immude[immudeod,])
totalname<- Immune[Immune$Gene_ID %in% rownames(total),]
totalname<- dplyr::select(totalname, "Gene_ID","Blastp_target_full_name")
total1<- cbind(total,Immune[match(rownames(total),Immune$Gene_ID),][,c(2,8)])

}


