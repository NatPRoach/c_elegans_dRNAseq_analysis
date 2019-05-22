#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
library(here)
theme_set(theme_cowplot())
setwd(here())

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}

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
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.justification = c(0.5, 0))
pdf("figures/figure4/figure4D_right.pdf",width=3,height = 2,colormodel="rgb")
p1
dev.off()
