#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(here)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}

setwd(here())

# df <- read.table("results/scratch/PASanalysis/L1_PAS_vs_polya.txt",header = TRUE,sep="\t")
# p <- ggplot(data=df,aes(x=PAS,y=length))+
#   geom_hline(yintercept=25,color="gray60")+
#   geom_hline(yintercept=50,color="gray60")+
#   geom_hline(yintercept=100,color="gray60")+
#   geom_violin(aes(fill=polya_color),color="black",show.legend = FALSE,draw_quantiles = c(0.5)) +
#   scale_y_continuous(expand=c(0,0))+
#   coord_cartesian(ylim=c(0,200))+
#   scale_color_manual(values=c("#a7cee2","#a7cee2","#a7cee2"))+
#   scale_fill_manual(values=c("#a7cee2","#a7cee2","#a7cee2"))+
#   xlab("PAS")+
#   ylab("polyA tail lengths")+
#   theme(text = element_text(size = 10,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=10,family = "Helvetica"),
#         axis.text.y = element_text(size=10,family = "Helvetica"),
#         legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
#         legend.position = "top",
#         legend.title = element_blank(),
#         legend.justification = c(0.5, 0))
# 
# 
# pdf(file="figures/figure4/figure4B.pdf",height=2.5,width=4.1667,colormodel="rgb")
# print(p)
# dev.off()
# 
# wilcox.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="altPAS",]$length)$p.value
# wilcox.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="noPAS",]$length)$p.value
# wilcox.test(df[df$PAS=="altPAS",]$length,df[df$PAS=="noPAS",]$length)$p.value
# ks.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="altPAS",]$length)$p.value
# ks.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="noPAS",]$length)$p.value
# ks.test(df[df$PAS=="altPAS",]$length,df[df$PAS=="noPAS",]$length)$p.value
# 
# median(df[df$PAS=="AATAAA",]$length)
# median(df[df$PAS=="altPAS",]$length)
# median(df[df$PAS=="noPAS",]$length)

### Sensitive

df <- read.table("results/scratch/PASanalysis/L1_sensitive_PAS_vs_polya.txt",header = TRUE,sep="\t")
p <- ggplot(data=df,aes(x=PAS,y=length))+
  geom_hline(yintercept=25,color="gray60")+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=100,color="gray60")+
  geom_violin(aes(fill=polya_color),color="black",show.legend = FALSE,draw_quantiles = c(0.5)) +
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,200))+
  scale_color_manual(values=c("#a7cee2","#a7cee2","#a7cee2"))+
  scale_fill_manual(values=c("#a7cee2","#a7cee2","#a7cee2"))+
  xlab("PAS")+
  ylab("polyA tail lengths")+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.justification = c(0.5, 0))


pdf(file="figures/figure4/figure4B_sensitive.pdf",height=2.5,width=4.1667,colormodel="rgb")
print(p)
dev.off()

wilcox.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="altPAS",]$length)$p.value
wilcox.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="noPAS",]$length)$p.value
wilcox.test(df[df$PAS=="altPAS",]$length,df[df$PAS=="noPAS",]$length)$p.value
ks.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="altPAS",]$length)$p.value
ks.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="noPAS",]$length)$p.value
ks.test(df[df$PAS=="altPAS",]$length,df[df$PAS=="noPAS",]$length)$p.value

median(df[df$PAS=="AATAAA",]$length)
median(df[df$PAS=="altPAS",]$length)
median(df[df$PAS=="noPAS",]$length)

### Stringent

df <- read.table("results/scratch/PASanalysis/L1_stringent_PAS_vs_polya.txt",header = TRUE,sep="\t")
p <- ggplot(data=df,aes(x=PAS,y=length))+
  geom_hline(yintercept=25,color="gray60")+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=100,color="gray60")+
  geom_violin(aes(fill=polya_color),color="black",show.legend = FALSE,draw_quantiles = c(0.5)) +
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,200))+
  scale_color_manual(values=c("#a7cee2","#a7cee2","#a7cee2"))+
  scale_fill_manual(values=c("#a7cee2","#a7cee2","#a7cee2"))+
  xlab("PAS")+
  ylab("polyA tail lengths")+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.justification = c(0.5, 0))


pdf(file="figures/figure4/figure4B_stringent.pdf",height=2.5,width=4.1667,colormodel="rgb")
print(p)
dev.off()

wilcox.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="altPAS",]$length)$p.value
wilcox.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="noPAS",]$length)$p.value
wilcox.test(df[df$PAS=="altPAS",]$length,df[df$PAS=="noPAS",]$length)$p.value
ks.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="altPAS",]$length)$p.value
ks.test(df[df$PAS=="AATAAA",]$length,df[df$PAS=="noPAS",]$length)$p.value
ks.test(df[df$PAS=="altPAS",]$length,df[df$PAS=="noPAS",]$length)$p.value

median(df[df$PAS=="AATAAA",]$length)
median(df[df$PAS=="altPAS",]$length)
median(df[df$PAS=="noPAS",]$length)
