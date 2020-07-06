#!/usr/bin/Rscript

library(tximport)
library(readr)
library(tidyr)
library(dplyr)
library(cowplot)
library(tibble)
library(DESeq2)

dir <- "./"

run <- c("cm_c_1_S1", "cm_c_2_S9", "cm_t_1_S2", "cm_t_2_S10")
files <- file.path(dir, run, "quant.sf")
names(files) <- run
all(file.exists(files))
tx2gene <- read_tsv("cm_tx2gene.tsv")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = c("species", "treatment"), sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- gsub("^c$", "control", sampleTable$treatment)
sampleTable$treatment <- gsub("^t$", "treatment", sampleTable$treatment)
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, -library, -species)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment)

bpgene<- read.csv("../cm_full_tax_nr_metazoa",header=F,sep="\t")
dds<- dds[rownames(dds) %in% bpgene$V1,]
dds<- dds[rowSums(counts(dds)>0) > 1, ]

dds <- DESeq(dds, parallel = TRUE, fitType = "local")
resLFC <- lfcShrink(dds, coef=2,type="apeglm")
png("./cm_MAplot.png")
plotMA(resLFC, ylim=c(-7,7))
dev.off()

rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
pdf("./pca.pdf")
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
pdf("./cm_PCA.pdf",21,7)
plot_grid(pc12, pc13, pc23, nrow = 1)
dev.off()

#res <- results(dds, alpha = 0.01) %>% subset(padj < 0.01 & abs(log2FoldChange) > 1)
res_tvsc <- results(dds, alpha = 0.01, tidy = TRUE, contrast = c("treatment", "treatment", "control"))
res_TvsC_sig_up <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.01 & log2FoldChange > 1)
write_csv(res_TvsC_sig_up,"./cmres_TvsC_sig_up.csv")
res_TvsC_sig_down <- dplyr::filter(as.data.frame(res_tvsc), padj < 0.01 & log2FoldChange < -1)
write_csv(res_TvsC_sig_down,"./cmres_TvsC_sig_down.csv")

#if (False){
Immune<-read.table("../../immune/cm_immune_curate_nrtaxcont",header=F,sep="\t")
colnames(Immune)<- c("clust","Gene_ID","Family","hmm1","hmm2","iso","type","len","val","protid","blaste","blastname","blastiden","fullname","tax","pfamd","pfamdname")
Immune_up<- res_TvsC_sig_up[res_TvsC_sig_up$row %in% Immune$Gene_ID,]
Immune_down<- res_TvsC_sig_down[res_TvsC_sig_down$row %in% Immune$Gene_ID,]
library(pheatmap)
immupe<- assay(rld)[as.character(Immune_up[,c("row")]),]
immude<- assay(rld)[as.character(Immune_down[,c("row")]),]
immupeod<- order(rowSums(immupe[,c(3,4)]),decreasing=T)
immudeod<- order(rowSums(immude[,c(3,4)]),decreasing=T)

refe<- assay(rld)["TRINITY_DN131018_c0_g1",]

total<- rbind(refe,immupe[immupeod,], immude[immudeod,])
#totalname<- Immune[Immune$Gene_ID %in% rownames(total),]
#totalname<- dplyr::select(totalname, "Gene_ID","Blastp_target_full_name")
#total1<- cbind(total,Immune[match(rownames(total),Immune$Gene_ID),][,c(2,8)])
pdf("./cm_regulate_immune_heatmap.pdf",10,14)
pheatmap(total,cluster_rows=F, annotation_col=sampleTable)
dev.off()
#}
