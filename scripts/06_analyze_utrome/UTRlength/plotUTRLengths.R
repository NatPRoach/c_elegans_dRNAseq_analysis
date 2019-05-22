#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
#library(plyr)
#library(ggpubr)
library(here)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  our_utrs=args[1]
}else{
  our_utrs='#DFC27D'
}
setwd(here())

df1 <- read.table("results/scratch/UTRlength/utr_lengths.txt",sep = "\t",header = TRUE )

g1 <- ggplot(df1,aes(x=reorder(stage,x_order),y=utr_length))+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=150,color="gray60")+
  geom_hline(yintercept=250,color="gray60")+
  geom_violin(draw_quantiles = c(0.50),color="black",fill=our_utrs)+
  ylab("3'UTR Length")+
  scale_x_discrete(labels=c("L1","L2","L3","L4","yAd", "mAd", "male","all"))+
  xlab("")+
  expand_limits(y=0)+
  scale_y_continuous(limits=c(0,500),expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family = "Helvetica"),
        legend.title = element_blank())

pdf("figures/figure3/figure3E.pdf",height=2.5,width=4.25,colormodel="rgb")
g1
dev.off()
# 
# ks.test(df1$utr_length[df1$stage=="L1"],df1$utr_length[df1$stage=="L2"])
# ks.test(df1$utr_length[df1$stage=="L2"],df1$utr_length[df1$stage=="L3"])
# ks.test(df1$utr_length[df1$stage=="L3"],df1$utr_length[df1$stage=="L4"])
# ks.test(df1$utr_length[df1$stage=="L4"],df1$utr_length[df1$stage=="young adult"])
# ks.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="mature adult"])
# ks.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="male"])
# 
# wilcox.test(df1$utr_length[df1$stage=="L1"],df1$utr_length[df1$stage=="L2"])
# wilcox.test(df1$utr_length[df1$stage=="L2"],df1$utr_length[df1$stage=="L3"])
# wilcox.test(df1$utr_length[df1$stage=="L3"],df1$utr_length[df1$stage=="L4"])
# wilcox.test(df1$utr_length[df1$stage=="L4"],df1$utr_length[df1$stage=="young adult"])
# wilcox.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="mature adult"])
# wilcox.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="male"])