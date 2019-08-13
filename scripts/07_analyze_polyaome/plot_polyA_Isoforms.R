#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(here)
setwd(here())

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}

# tmp <- read.table("results/scratch/polya/fulls/all.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# All <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
#
# p1<- ggplot(All) +
#   geom_hline(yintercept=25,color="gray60")+
#   geom_hline(yintercept=50,color="gray60")+
#   geom_hline(yintercept=100,color="gray60")+
#   geom_violin(aes(x="Y37E3.8b.1",y=Y37E3.8b.1,fill=""),draw_quantiles = 0.5) +
#   geom_violin(aes(x="Y37E3.8a",y=Y37E3.8a,fill=""),draw_quantiles = 0.5) +
#   scale_fill_manual(values=c(polya_color,polya_color))+
#   scale_y_continuous(expand=c(0,0))+
#   ylab("estimated poly(A) length")+
#   xlab("")+
#   theme(text = element_text(size = 10,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=10,family = "Helvetica"),
#         axis.text.y = element_text(size=10,family = "Helvetica"),
#         legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.justification = c(0.5, 0))
# pdf("figures/figure4/figure4C_right.pdf",width=3,height = 2,colormodel="rgb")
# p1
# dev.off()

### Stringent

tmp <- read.table("results/scratch/polya/fulls/all_stringent.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
All <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
length(All[!is.na(All$Y37E3.8b.1),]$Y37E3.8b.1)
length(All[!is.na(All$Y37E3.8a),]$Y37E3.8a)

p1<- ggplot(All) +
  # geom_hline(yintercept=25,color="gray60")+
  # geom_hline(yintercept=50,color="gray60")+
  # geom_hline(yintercept=100,color="gray60")+
  geom_hline(yintercept=25,colour="#990000", linetype="dashed")+
  geom_hline(yintercept=50,colour="#990000", linetype="dashed")+
  geom_hline(yintercept=100,colour="#990000", linetype="dashed")+
  geom_violin(aes(x="Y37E3.8b.1",y=Y37E3.8b.1,fill=""),draw_quantiles = 0.5) +
  geom_violin(aes(x="Y37E3.8a",y=Y37E3.8a,fill=""),draw_quantiles = 0.5) +
  scale_fill_manual(values=c(polya_color,polya_color))+
  scale_y_continuous(expand=c(0,0))+
  ylab("estimated poly(A) length")+
  xlab("")+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.justification = c(0.5, 0))
# pdf("figures/figure4/figure4C_right_stringent.pdf",width=3,height = 2,colormodel="rgb")
pdf("figures/figure4/figure4C_right.pdf",width=3,height = 2,colormodel="rgb")
p1
dev.off()

# ### Sensitive
# 
# tmp <- read.table("results/scratch/polya/fulls/all_sensitive.full.polya",header = FALSE,fill=TRUE,stringsAsFactors=FALSE)
# All <- setNames(data.frame(t(tmp[,-1])),tmp[,1])
# 
# p1<- ggplot(All) +
#   # geom_hline(yintercept=25,color="gray60")+
#   # geom_hline(yintercept=50,color="gray60")+
#   # geom_hline(yintercept=100,color="gray60")+
#   geom_hline(yintercept=25,colour="#990000", linetype="dashed")+
#   geom_hline(yintercept=50,colour="#990000", linetype="dashed")+
#   geom_hline(yintercept=100,colour="#990000", linetype="dashed")+
#   geom_violin(aes(x="Y37E3.8b.1",y=Y37E3.8b.1,fill=""),draw_quantiles = 0.5) +
#   geom_violin(aes(x="Y37E3.8a",y=Y37E3.8a,fill=""),draw_quantiles = 0.5) +
#   scale_fill_manual(values=c(polya_color,polya_color))+
#   scale_y_continuous(expand=c(0,0))+
#   ylab("estimated poly(A) length")+
#   xlab("")+
#   theme(text = element_text(size = 10,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=10,family = "Helvetica"),
#         axis.text.y = element_text(size=10,family = "Helvetica"),
#         legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
#         legend.position = "none",
#         legend.title = element_blank(),
#         legend.justification = c(0.5, 0))
# pdf("figures/figure4/figure4C_right_sensitive.pdf",width=3,height = 2,colormodel="rgb")
# p1
# dev.off()