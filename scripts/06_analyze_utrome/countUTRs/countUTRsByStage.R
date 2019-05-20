#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
library(here)
theme_set(theme_cowplot())

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  mangone_utrs=args[1]
  jan_utrs=args[2]
  our_utrs=args[3]
  genes_with_novel_utrs=args[4]
  novel_utrs=args[5]
}else{
  mangone_utrs='#B3DE8E'
  jan_utrs='#F99B9B'
  our_utrs='#A7CEE2'
  genes_with_novel_utrs="#2579B2"
  novel_utrs="#389E34"
}

setwd(here())
df <- read.table("results/scratch/countUTRs/utr_counts.txt",header=TRUE,sep='\t')
# df$dataset <- factor(df$dataset,levels = c("Our genes","Waterston genes","Our isoforms","Waterston isoforms"))
g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(jan_utrs,mangone_utrs,our_utrs),labels=c("Our UTRs" = "This\nstudy","Mangone UTRs"="Mangone\n3'UTRs","Jan UTRs"="Jan\n3'UTRs"))+
  ylab("Number 3' UTRs identified")+
  xlab("Stage")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
        #legend.key.size = unit(1.25, 'lines'),
        #legend.position = "top",
        legend.position = c(0.4,0.85),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        #legend.justification = c(0.5, 0),
        legend.justification = c(0.5, 0))
pdf(file="figures/figure3/figure3A.pdf",height=2.5,width=4.25,colormodel="rgb")
print(g2)
dev.off()

# setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/novelUTRs/")
df <- read.table("results/scratch/novelUTRs/novelUTRs.txt",header=TRUE,sep='\t')
# df$dataset <- factor(df$dataset,levels = c("Our genes","Waterston genes","Our isoforms","Waterston isoforms"))
g3 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(genes_with_novel_utrs,novel_utrs),labels=c("genes with novel UTRs"="genes with\nnovel 3'UTRs","novel UTRs" = "novel 3'UTRs"))+
  ylab("Number identified")+
  xlab("Stage")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
        #legend.key.size = unit(1.25, 'lines'),
        #legend.position = "top",
        legend.position = c(0.5,0.85),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        #legend.justification = c(0.5, 0),
        legend.justification = c(0.5, 0))
pdf(file="figures/figure3/figure3C.pdf",height=2.5,width=4.25,colormodel="rgb")
print(g3)
dev.off()
