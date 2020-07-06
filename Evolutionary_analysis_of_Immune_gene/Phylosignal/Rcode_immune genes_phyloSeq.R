#!/usr/bin/Rscript

library(lattice)
library(Rmisc)
library(ggplot2)
library(readxl)
termite <- read_excel("/media/shulinhe/DATA/Full_analysis/Phylosignal/termite_phyloseq_data.xlsx")
View(termite)
term <- summarySE(termite, measurevar="genes", groupvars=c("sociality"))

as.numeric(termite$busco)
as.factor(termite$sociality)

plot(termite$busco, termite$genes)
cor.test(termite$genes, termite$busco)

level_order <- c('sol', 'sub', 'soc')

ggplot(termite, aes(x=factor(sociality, level = level_order), y=genes, colour=sociality)) +
  geom_point(colour="white", shape=21, size = 4, 
  aes(fill = factor(sociality))) + 
  scale_fill_manual(values=c("blue", "cyan4", "green"))

ggplot(term, aes(x=factor(sociality, level = level_order), y=genes, colour=sociality)) + 
  geom_errorbar(aes(ymin=genes-se, ymax=genes+se), width=.1) +
  geom_line() +
  geom_point()

ggplot(termite, aes(x=busco, y=genes, color=sociality)) + 
  geom_point(colour="white", shape=21, size = 4, 
             aes(fill = factor(sociality))) + 
  scale_fill_manual(values=c("blue", "cyan4", "green"))

ggplot(termite, aes(x=busco, y=genes, color=sociality)) + 
  geom_point(shape=1) +
  scale_colour_hue(l=50) +
  geom_smooth(method=lm, se=FALSE)

###TEST PHYLOGENETIC COMPONENT WITH PHYLOSIGNAL###

library(phylosignal)
library(phylobase)
library(ape)

#ultrametric
text.string<-
  "(Bg:13,(Bo:12,((Cm:1,Ca:1):10,(Md:10,(Zn:9,((Kf:2,(Nc:1,Cb:1):1):6,(Pi:7,((Cf:1,Rf:1):5,(Ms:5,(LS17:4,(LS29:3,((LS19:1,PG24:1):1,LS26:2):1):1):1):1):1):1):1):1):1):1):1);"
tre<-read.tree(text=text.string)
plot(tre)
dat <- list()
dat$busco <- termite$busco
dat$genes <- termite$genes
dat <- as.data.frame(dat)
termp4d <- phylo4d(tre, dat)
barplot.phylo4d(termp4d, tree.type = "phylo", tree.ladderize = TRUE)
phyloSignal(p4d = termp4d, method = "all")

busco.cg <- phyloCorrelogram(termp4d, trait = "busco")
genes.cg <- phyloCorrelogram(termp4d, trait = "genes")
plot(busco.cg)
plot(genes.cg)

#ladder
text.string1<-
  "(Bg:1,(Bo:1,((Cm:1,Ca:1):1,(Md:1,(Zn:1,((Kf:1,(Nc:1,Cb:1):1):1,(Pi:1,((Cf:1,Rf:1):1,(Ms:1,(LS17:1,(LS29:1,((LS19:1,PG24:1):1,LS26:1):1):1):1):1):1):1):1):1):1):1):1);"
tre1<-read.tree(text=text.string1)
plot(tre1)
dat1 <- list()
dat1$busco <- termite$busco
dat1$genes <- termite$genes
dat1 <- as.data.frame(dat1)
termp4d1 <- phylo4d(tre1, dat1)
barplot.phylo4d(termp4d1, tree.type = "phylo", tree.ladderize = TRUE)
phyloSignal(p4d = termp4d1, method = "all")

busco.cg1 <- phyloCorrelogram(termp4d1, trait = "busco")
genes.cg1 <- phyloCorrelogram(termp4d1, trait = "genes")
plot(busco.cg1)
plot(genes.cg1)

#actual tree
text.string2<-
  "(Bg:271.8456,(Bo:216.6575,((Ca:13.4524,Cm:13.4524):165.984,(Md:155.341,(Zn:147.4815,((Kf:67.2899,(Nc:43.943,Cb:43.943):23.3469):65.9293,(Pi:86.8775,((Cf:58.836,Rf:58.836):17.7823,(Ms:58.9308,(LS17:49.4615,(LS29:38.2166,((PG24:30.8428,LS19:30.8428):3.2148,LS26:34.0576):4.1589):11.2448):9.4693):17.6874):10.2592):46.3417):14.2623):7.8596):24.0949):37.2212):55.188);"

tre2<-read.tree(text=text.string2)
plot(tre2)
dat2 <- list()
dat2$busco <- termite$busco
dat2$genes <- termite$genes
dat2 <- as.data.frame(dat2)
termp4d2 <- phylo4d(tre2, dat2)
barplot.phylo4d(termp4d2, tree.type = "phylo", tree.ladderize = TRUE)
phyloSignal(p4d = termp4d2, method = "all")
###
#$stat
#           Cmean           I         K    K.star       Lambda
#busco 0.05849514 -0.05907386 0.3707184 0.4890398 0.0000695486
#genes 0.44913515  0.05502950 1.3905187 0.8689007 0.8302625689
#
#$pvalue
#      Cmean     I     K K.star      Lambda
#busco 0.178 0.467 0.365  0.286 1.000000000
#genes 0.002 0.023 0.002  0.008 0.007956169
###

busco.cg2 <- phyloCorrelogram(termp4d2, trait = "busco")
genes.cg2 <- phyloCorrelogram(termp4d2, trait = "genes")
plot(busco.cg2)
plot(genes.cg2)
