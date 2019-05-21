#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
#library(ggridges)
#library(lubridate)
#library(colorspace)
#library("dviz.supp")
library(gdata)
# library("cowplot")
# library("qvalue")
library(grid)
library(gridExtra)
library(here)
# setwd("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/")
setwd(here())
# median.quartile <- function(x){
#   out <- quantile(x, probs = c(0.25,0.5,0.75))
#   names(out) <- c("ymin","y","ymax")
#   return(out) 
# }

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}

##############################
#### Initialize color schemes for figures
##############################
# COL <- (c(rgb(187,240,255,maxColorValue=255),rgb(119,226,255,maxColorValue=255),rgb(0,132,169,maxColorValue=255),rgb(0,99,127,maxColorValue=255),rgb(0,66,85,maxColorValue=255)))
# COL2 <- c("#FF050D","#FF050D7F","#8587FF","#8587FF7F","#00C52D50","#00C52D7F")
# COL_dRNA_cDNA <- c("#AAC78E","#AAC78E7F","#1F78B4","#1F78B47F")

# tmp <- read.table("fulls/l1.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# L1 <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
# tmp <- read.table("fulls/l2.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# L2 <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
# tmp <- read.table("fulls/l3.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# L3 <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
# tmp <- read.table("fulls/l4.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# L4 <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
# tmp <- read.table("fulls/ya.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# YA <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
# tmp <- read.table("fulls/ga.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# GA <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
# tmp <- read.table("fulls/ml.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# ML <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
tmp <- read.table("results/scratch/polya/fulls/all.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
All <- setNames(data.frame(t(tmp[,-1])),tmp[,1])

# L1_sig <- read.table("fulls/l1.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# L1_comparison_count <- nrow(L1_sig)
# tmp <- qvalue(L1_sig$p.val)
# L1_sig$q.val <- tmp$qvalues
# L1_sig <- L1_sig[L1_sig$q.val < 0.05,]
# L1_sig_count <- nrow(L1_sig)
# L1_percent <- 100* L1_sig_count / L1_comparison_count
# 
# L2_sig <- read.table("fulls/l2.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# L2_comparison_count <- nrow(L2_sig)
# tmp <- qvalue(L2_sig$p.val)
# L2_sig$q.val <- tmp$qvalues
# L2_sig <- L2_sig[L2_sig$q.val < 0.05,]
# L2_sig_count <- nrow(L2_sig)
# L2_percent <- 100* L2_sig_count / L2_comparison_count
# 
# 
# L3_sig <- read.table("fulls/l3.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# L3_comparison_count <- nrow(L3_sig)
# tmp <- qvalue(L3_sig$p.val)
# L3_sig$q.val <- tmp$qvalues
# L3_sig <- L3_sig[L3_sig$q.val < 0.05,]
# L3_sig_count <- nrow(L3_sig)
# L3_percent <- 100* L3_sig_count / L3_comparison_count
# 
# 
# L4_sig <- read.table("fulls/l4.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# L4_comparison_count <- nrow(L4_sig)
# tmp <- qvalue(L4_sig$p.val)
# L4_sig$q.val <- tmp$qvalues
# L4_sig <- L4_sig[L4_sig$q.val < 0.05,]
# L4_sig_count <- nrow(L4_sig)
# L4_percent <- 100* L4_sig_count / L4_comparison_count
# 
# 
# YA_sig <- read.table("fulls/ya.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# YA_comparison_count <- nrow(YA_sig)
# tmp <- qvalue(YA_sig$p.val)
# YA_sig$q.val <- tmp$qvalues
# YA_sig <- YA_sig[YA_sig$q.val < 0.05,]
# YA_sig_count <- nrow(YA_sig)
# YA_percent <- 100* YA_sig_count / YA_comparison_count
# 
# 
# GA_sig <- read.table("fulls/ga.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# GA_comparison_count <- nrow(GA_sig)
# tmp <- qvalue(GA_sig$p.val)
# GA_sig$q.val <- tmp$qvalues
# GA_sig <- GA_sig[GA_sig$q.val < 0.05,]
# GA_sig_count <- nrow(GA_sig)
# GA_percent <- 100* GA_sig_count / GA_comparison_count
# 
# 
# ML_sig <- read.table("fulls/ml.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# ML_comparison_count <- nrow(ML_sig)
# tmp <- qvalue(ML_sig$p.val)
# ML_sig$q.val <- tmp$qvalues
# ML_sig <- ML_sig[ML_sig$q.val < 0.05,]
# ML_sig_count <- nrow(ML_sig)
# ML_percent <- 100* ML_sig_count / ML_comparison_count
# 
# 
# All_sig <- read.table("fulls/all.sig.polya",header=TRUE,stringsAsFactors=FALSE)
# All_comparison_count <- nrow(All_sig)
# tmp <- qvalue(All_sig$p.val)
# All_sig$q.val <- tmp$qvalues
# All_sig <- All_sig[All_sig$q.val < 0.05,]
# All_sig_count <- nrow(All_sig)
# All_percent <- 100* All_sig_count / All_comparison_count



# p1<- ggplot(L1) +
#   geom_violin(aes(x="W09C5.6a",y=W09C5.6a,fill=COL2[3])) +
#   geom_violin(aes(x="W09C5.6b.5",y=W09C5.6b.5,fill=COL2[4])) +
#   geom_jitter(aes(x="W09C5.6a",y=W09C5.6a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="W09C5.6b.5",y=W09C5.6b.5),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("W09C5.6a","W09C5.6b.5"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   ggtitle(expression(italic("rpl-31")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("polya_vs_isoform2.pdf",plot = p1,width=15,height = 10,units = "in")
# 
# p2<- ggplot(GA) +
#   geom_violin(aes(x="F28C6.7a",y=F28C6.7a,fill=COL2[3])) +
#   geom_violin(aes(x="F28C6.7b",y=F28C6.7b,fill=COL2[4])) +
#   geom_jitter(aes(x="F28C6.7a",y=F28C6.7a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="F28C6.7b",y=F28C6.7b),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("F28C6.7a","F28C6.7b"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   ggtitle(expression(italic("rpl-26")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("polya_vs_isoform3.pdf",plot = p2,width=15,height = 10,units = "in")
# 
# 
# p3<- ggplot(YA) +
#   geom_violin(aes(x="Y37E3.8a",y=Y37E3.8a,fill=COL2[3])) +
#   geom_violin(aes(x="Y37E3.8b.1",y=Y37E3.8b.1,fill=COL2[4])) +
#   geom_jitter(aes(x="Y37E3.8a",y=Y37E3.8a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="Y37E3.8b.1",y=Y37E3.8b.1),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("Y37E3.8a","Y37E3.8b.1"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   ggtitle("WBGene00021350")+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("polya_vs_isoform4.pdf",plot = p3,width=15,height = 10,units = "in")
# 
# # M117.2a.2
# # M117.2a.3
# 
# p4<- ggplot(L3) +
#   geom_violin(aes(x="M117.2a.2",y=M117.2a.2,fill=COL2[3])) +
#   geom_violin(aes(x="M117.2a.3",y=M117.2a.3,fill=COL2[4])) +
#   geom_jitter(aes(x="M117.2a.2",y=M117.2a.2),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="M117.2a.3",y=M117.2a.3),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("M117.2a.2","M117.2a.3"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   ggtitle(expression(italic("par-5")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("polya_vs_isoform4.pdf",plot = p4,width=15,height = 10,units = "in")

p1<- ggplot(All) +
  geom_hline(yintercept=25,color="gray60")+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=100,color="gray60")+
  geom_violin(aes(x="Y37E3.8b.1",y=Y37E3.8b.1,fill=""),draw_quantiles = 0.5) +
  geom_violin(aes(x="Y37E3.8a",y=Y37E3.8a,fill=""),draw_quantiles = 0.5) +
  #geom_jitter(aes(x="Y37E3.8b.1",y=Y37E3.8b.1),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
  #geom_jitter(aes(x="Y37E3.8a",y=Y37E3.8a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
  #scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("W09C5.6a","W09C5.6b.5"))+
  scale_fill_manual(values=c(polya_color,polya_color))+
  scale_y_continuous(expand=c(0,0))+
  ylab("estimated poly(A) length")+
  xlab("")+
  #ggtitle(expression(italic("rpl-31")))+
  #ylim(c(0,300))+
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
pdf("figures/figure4/figure4C_right.pdf",width=3,height = 2,colormodel="rgb")
p1
dev.off()
# ggsave("plots/polya_vs_isoform.pdf",plot = p1,width=15,height = 10,units = "in")


# p1<- ggplot(L1) +
#   geom_violin(aes(x="W09C5.6a",y=W09C5.6a),fill="#A7CEE2",color="black",draw_quantiles = 0.5) +
#   geom_violin(aes(x="W09C5.6b.5",y=W09C5.6b.5),fill="#A7CEE2",color="black",draw_quantiles = 0.5) +
#   #geom_jitter(aes(x="W09C5.6a",y=W09C5.6a),alpha = 0.01,shape=16, position=position_jitter(0.2)) +
#   #geom_jitter(aes(x="W09C5.6b.5",y=W09C5.6b.5),alpha = 0.01,shape=16, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",labels=c("W09C5.6a","W09C5.6b.5"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   #ggtitle(expression(italic("rpl-31")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size = 8,family="Helvetica"),
#         axis.text.x = element_text(size=8,family="Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size=8,family="Helvetica"))
# 
# pdf("plots/polya_vs_isoform.pdf",width=3.5,height = 2.5,colormodel="rgb")
# p1
# dev.off()
# 
# 
# p1<- ggplot(L1) +
#   geom_violin(aes(x="Y37E3.8b.1",y=Y37E3.8b.1,fill=COL2[3])) +
#   geom_violin(aes(x="Y37E3.8a",y=Y37E3.8a,fill=COL2[4])) +
#   geom_jitter(aes(x="Y37E3.8b.1",y=Y37E3.8b.1),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="Y37E3.8a",y=Y37E3.8a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("W09C5.6a","W09C5.6b.5"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   #ggtitle(expression(italic("rpl-31")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),plot.title = element_text(size=26),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("plots/polya_vs_isoform.pdf",plot = p1,width=15,height = 10,units = "in")
# 
# 
# p1<- ggplot(L1) +
#   geom_violin(aes(x="W09C5.6a",y=W09C5.6a,fill=COL2[3])) +
#   geom_violin(aes(x="W09C5.6b.5",y=W09C5.6b.5,fill=COL2[4])) +
#   geom_jitter(aes(x="W09C5.6a",y=W09C5.6a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="W09C5.6b.5",y=W09C5.6b.5),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("W09C5.6a","W09C5.6b.5"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   ggtitle(expression(italic("rpl-31")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),plot.title = element_text(size=26),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("plots/polya_vs_isoform.pdf",plot = p1,width=15,height = 10,units = "in")
# 
# p2<- ggplot(L1) +
#   geom_violin(aes(x="Y71F9AL.13",y=Y71F9AL.13,fill=COL2[3])) +
#   geom_violin(aes(x="Y71F9AL.13a",y=Y71F9AL.13a,fill=COL2[4])) +
#   geom_violin(aes(x="Y71F9AL.13b.7",y=Y71F9AL.13b.7,fill=COL2[5])) +
#   geom_jitter(aes(x="Y71F9AL.13",y=Y71F9AL.13),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="Y71F9AL.13a",y=Y71F9AL.13a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="Y71F9AL.13b.7",y=Y71F9AL.13b.7),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:5],labels=c("Y71F9AL.13","Y71F9AL.13a","Y71F9AL.13b.7"))+
#   ylab("")+
#   xlab("")+
#   ggtitle(expression(italic("rpl-1")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),plot.title = element_text(size=26),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("plots/polya_vs_isoform2.pdf",plot = p2,width=15,height = 10,units = "in")
# 
# 
# p3<- ggplot(L1) +
#   geom_violin(aes(x="JC8.3a",y=JC8.3a,fill=COL2[3])) +
#   geom_violin(aes(x="JC8.3",y=JC8.3,fill=COL2[4])) +
#   geom_violin(aes(x="JC8.3c.4",y=JC8.3c.4,fill=COL2[5])) +
#   geom_jitter(aes(x="JC8.3a",y=JC8.3a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="JC8.3",y=JC8.3),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="JC8.3c.4",y=JC8.3c.4),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:5],labels=c("JC8.3a","JC8.3","JC8.3c.4"))+
#   ylab("")+
#   xlab("")+
#   ggtitle(expression(italic("rpl-12")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),plot.title = element_text(size=26),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("plots/polya_vs_isoform3.pdf",plot = p3,width=15,height = 10,units = "in")
# 
# p4<- ggplot(L1) +
#   geom_violin(aes(x="M18.7a",y=M18.7a,fill=COL2[3])) +
#   geom_violin(aes(x="M18.7b",y=M18.7b,fill=COL2[4])) +
#   geom_jitter(aes(x="M18.7a",y=M18.7a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="M18.7b",y=M18.7b),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("M18.7a","M18.7b"))+
#   ylab("Estimated Poly(A) Length")+
#   xlab("")+
#   ggtitle(expression(italic("aly-3")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),plot.title = element_text(size=26),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("plots/polya_vs_isoform4.pdf",plot = p4,width=15,height = 10,units = "in")
# 
# p5<- ggplot(L1) +
#   geom_violin(aes(x="Y37E3.8a",y=Y37E3.8a,fill=COL2[3])) +
#   geom_violin(aes(x="Y37E3.8b.1",y=Y37E3.8b.1,fill=COL2[4])) +
#   geom_violin(aes(x="Y37E3.8b.2",y=Y37E3.8b.2,fill=COL2[5])) +
#   geom_jitter(aes(x="Y37E3.8a",y=Y37E3.8a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="Y37E3.8b.1",y=Y37E3.8b.1),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="Y37E3.8b.2",y=Y37E3.8b.2),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:5],labels=c("Y37E3.8a","Y37E3.8b.1","Y37E3.8b.2"))+
#   ylab("")+
#   xlab("")+
#   ggtitle(expression(italic("Y37E3.8")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),plot.title = element_text(size=26),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("plots/polya_vs_isoform5.pdf",plot = p5,width=15,height = 10,units = "in")
# 
# p6<- ggplot(L1) +
#   geom_violin(aes(x="R07B7.3a",y=R07B7.3a,fill=COL2[3])) +
#   geom_violin(aes(x="R07B7.3b",y=R07B7.3b,fill=COL2[4])) +
#   geom_jitter(aes(x="R07B7.3a",y=R07B7.3a),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   geom_jitter(aes(x="R07B7.3b",y=R07B7.3b),alpha = 0.2,shape=16,size=2.5, position=position_jitter(0.2)) +
#   scale_fill_discrete(name="Isoform",breaks=COL2[3:4],labels=c("R07B7.3a","R07B7.3b"))+
#   ylab("")+
#   xlab("")+
#   ggtitle(expression(italic("pqn-53")))+
#   ylim(c(0,300))+
#   theme(text = element_text(size=32),plot.title = element_text(size=26),axis.text.x = element_text(size=24),axis.text.y = element_text(size=28),legend.position="none")
# ggsave("plots/polya_vs_isoform6.pdf",plot = p6,width=15,height = 10,units = "in")
# 
# 
# 
# pdf(file="plots/six_examples.pdf",height=12,width=26,colormodel="rgb")
# grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2,ncol=3)
# dev.off()
# 
