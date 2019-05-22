#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(here)
library(scales)
args = commandArgs(trailingOnly=TRUE)



if (length(args)!=0){
splice_isoform_color = args[1]
utr_color = args[2]
fl_isoform_color = args[3]
}else{
  splice_isoform_color= "#C9B3D5"
  utr_color = "#67A9CF"
  fl_isoform_color = "#02818A"
}

setwd(here())


df <- read.table("results/scratch/countCombinatorialIsoforms/combinatorial_isoform_count.txt",sep = "\t", header= TRUE)

df$dataset <- factor(df$dataset,levels = c("splice isoforms","utrs","full length isoforms"))
g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
  geom_col(position=position_dodge(0.75),width=0.75,aes(fill = dataset),color="black")+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(splice_isoform_color,utr_color,fl_isoform_color),labels = c("splice isoforms"="splice\nisoforms","utrs" = "3'UTRs","full length isoforms"="full-length\nisoforms"))+
  ylab(bquote('Number identified (x'~ 10^3~')' ))+
  xlab("Stage")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0),labels=function(x)x/1000)+
  theme(text = element_text(size = 10,family="Helvetica"),
        axis.text.x = element_text(size=10,family="Helvetica"),
        axis.text.y = element_text(size=10,family="Helvetica"),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = c(0.5,0.80),
        legend.justification = c(0.5, 0),
        legend.direction = "horizontal",
        axis.ticks.x = element_blank())

pdf(file="figures/figure1/figure1F.pdf",height=2.0834,width=4.35,colormodel="rgb")
print(g2)
dev.off()


df2 <- read.table("results/scratch/countCombinatorialIsoforms/combinatorial_isoform_count_subsampled.txt",sep="\t",header=TRUE)
df2$stage <- factor(df2$stage,levels = c("L1","L2","L3","L4","young adult","mature adult","male","all"))


df3 <- df2[df2$stage != "all",]
g4 <- ggplot(df3,aes(x=read_count,y=counts,color=stage))+
  geom_line()+
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#a6761d","#b3b3b3","#000000"))+
  ylab("Number of full-length isoforms")+
  xlab("Number of reads retained")+
  expand_limits(y=0)+
  scale_x_continuous(labels=comma)+
  scale_y_continuous(limits =c(0,17000),expand=c(0,0))+
  theme(text = element_text(size = 10,family="Helvetica"),
        axis.text.x = element_text(size=10,family="Helvetica"),
        axis.text.y = element_text(size=10,family="Helvetica"),
        legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = c(0.5,0.1),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.justification = c(0.5, 0))

pdf(file="figures/figure1/figure1G.pdf",height=2.5,width=5,colormodel="rgb")
print(g4)
dev.off()

df3 <- df2[df2$stage == "all",]
g4 <- ggplot(df3,aes(x=read_count,y=counts,color=stage))+
  geom_line()+
  scale_color_manual(values=c("#000000","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#000000"))+
  ylab("Number of full-length isoforms")+
  xlab("Number of reads retained")+
  expand_limits(y=c(0,27000))+
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family="Helvetica"),
        axis.text.x = element_text(size=10,family="Helvetica"),
        axis.text.y = element_text(size=10,family="Helvetica"),
        legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = c(0.5,0.1),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.justification = c(0.5, 0))

pdf(file="figures/supplementals/sfigure2/sfigure2A.pdf",height=2.5,width=5,colormodel="rgb")
print(g4)
dev.off()

