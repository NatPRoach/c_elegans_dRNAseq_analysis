#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
library(here)
theme_set(theme_cowplot())

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}

# setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/")
setwd(here())
data <- read.table("results/scratch/polya/splice_v_intron.txt",header=FALSE,sep="\t")
colnames(data) <- c("estimated polyA length", "type")

p<- ggplot(data,aes(x=type,y=`estimated polyA length`))+
  geom_hline(yintercept=25,color="gray60")+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=100,color="gray60")+
  geom_violin(aes(fill=type),show.legend = FALSE,draw_quantiles = c(0.5))+
  scale_fill_manual(values=c(polya_color,polya_color))+
  # scale_y_continuous(expand=c(0,0),limits=c(0,300))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,300))+
  xlab("")+
  ylab("estimated poly(A) length")+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
        #legend.key.size = unit(1.25, 'lines'),
        legend.position = "top",
        legend.title = element_blank(),
        #legend.justification = c(0.5, 0),
        legend.justification = c(0.5, 0))
pdf(file="figures/figure4/figure4E.pdf",height=2.5,width=4.1667,colormodel="rgb")
p
dev.off()
