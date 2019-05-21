#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(qvalue)
library(ggplot2)
library(cowplot)
library(here)
theme_set(theme_cowplot())
# setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/polya_correlations/utr/")
setwd(here())

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}

# df2 <- read.table("all.utr_polya_sig.txt",sep="\t")
# colnames(df2) <- c("gene","utrs","pval")
# qvals2 <- qvalue(df2$pval)
# df2$qval <- qvals2$qvalues
# all_sig <- df2[df2$qval < 0.05,]
# 
# write.table(sig2$gene,file="utrVsTxGeneList.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 
# df3 <- read.table("l1.utr_polya_sig.txt",sep="\t")
# colnames(df3) <- c("gene","utrs","pval")
# qvals3 <- qvalue(df3$pval)
# df3$qval <- qvals3$qvalues
# l1_sig <- df3[df3$qval < 0.05,]

# l1 <- read.table("l1.utr_polya_lengths.txt",sep="\t",header = TRUE)
# df5 <- l1[l1$gene == "WBGene00006924",]
# 
# g1 <- ggplot(df5,aes(x=utr,y=length))+
#   geom_violin(fill="#A7CEE2",color="black",draw_quantiles = 0.5)+
#   scale_x_discrete(labels=c("WBGene00006924-cluster0"="vig-1 3'UTR 0","WBGene00006924-cluster1"="vig-1 3'UTR 1"))+
#   #geom_jitter(alpha=0.1)+
#   ylab("Estimated Poly(A) Length")+
#   theme(text = element_text(size = 8,family="Helvetica"),
#         axis.text.x = element_text(size=8,family="Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size=8,family="Helvetica"))
# 
# pdf("plots/polya_vs_UTR.pdf",width=3.5,height = 2.5,colormodel="rgb")
# g1
# dev.off()

all <- read.table("results/scratch/polya/all.utr_polya_lengths.txt",sep="\t",header = TRUE)
df5 <- all[all$gene == "WBGene00003920",]
df5 <- df5[df5$utr == "WBGene00003920-cluster2" | df5$utr == "WBGene00003920-cluster0",]
ks.test(df5[df5$utr == "WBGene00003920-cluster2",]$length,df5[df5$utr == "WBGene00003920-cluster0",]$length)
p1<- ggplot(df5,aes(x=utr,y=length,fill="")) +
  geom_hline(yintercept=25,color="gray60")+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=100,color="gray60")+
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_manual(values=c(polya_color,polya_color))+
  ylab("estimated poly(A) length")+
  xlab("")+
  scale_x_discrete(labels = c("WBGene00003920-cluster0" = "par-5 UTR0", "WBGene00003920-cluster2" = "par-5 UTR2"))+
  #ggtitle(expression(italic("rpl-31")))+
  #ylim(c(0,300))+
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
        #legend.key.size = unit(1.25, 'lines'),
        legend.position = "none",
        legend.title = element_blank(),
        #legend.justification = c(0.5, 0),
        legend.justification = c(0.5, 0))
pdf("figures/figure4/figure4D_right.pdf",width=3,height = 2,colormodel="rgb")
p1
dev.off()

# g1 <- ggplot(df5,aes(x=utr,y=length))+
#   geom_violin()+
#   geom_jitter(alpha=0.1)
# 
# df3 <- read.table("l2.utr_polya_sig.txt",sep="\t")
# colnames(df3) <- c("gene","utrs","pval")
# qvals3 <- qvalue(df3$pval)
# df3$qval <- qvals3$qvalues
# l2_sig <- df3[df3$qval < 0.05,]
# 
# df3 <- read.table("l3.utr_polya_sig.txt",sep="\t")
# colnames(df3) <- c("gene","utrs","pval")
# qvals3 <- qvalue(df3$pval)
# df3$qval <- qvals3$qvalues
# l3_sig <- df3[df3$qval < 0.05,]
# 
# df3 <- read.table("l4.utr_polya_sig.txt",sep="\t")
# colnames(df3) <- c("gene","utrs","pval")
# qvals3 <- qvalue(df3$pval)
# df3$qval <- qvals3$qvalues
# l4_sig <- df3[df3$qval < 0.05,]
# 
# df3 <- read.table("ya.utr_polya_sig.txt",sep="\t")
# colnames(df3) <- c("gene","utrs","pval")
# qvals3 <- qvalue(df3$pval)
# df3$qval <- qvals3$qvalues
# ya_sig <- df3[df3$qval < 0.05,]
# 
# df3 <- read.table("ga.utr_polya_sig.txt",sep="\t")
# colnames(df3) <- c("gene","utrs","pval")
# qvals3 <- qvalue(df3$pval)
# df3$qval <- qvals3$qvalues
# ga_sig <- df3[df3$qval < 0.05,]
# 
# df3 <- read.table("ml.utr_polya_sig.txt",sep="\t")
# colnames(df3) <- c("gene","utrs","pval")
# qvals3 <- qvalue(df3$pval)
# df3$qval <- qvals3$qvalues
# ml_sig <- df3[df3$qval < 0.05,]