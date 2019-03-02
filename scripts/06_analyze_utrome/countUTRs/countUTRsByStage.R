#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
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

setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/countUTRs/")
df <- read.table("utr_counts.txt",header=TRUE,sep='\t')
# df$dataset <- factor(df$dataset,levels = c("Our genes","Waterston genes","Our isoforms","Waterston isoforms"))
g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(jan_utrs,mangone_utrs,our_utrs),labels=c("Our UTRs" = "This\nstudy","Mangone UTRs"="Mangone\nUTRs","Jan UTRs"="Jan\nUTRs"))+
  ylab("Number 3' UTRs identified")+
  xlab("Stage")+
  theme(text = element_text(size = 8,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,family = "Helvetica"),
        axis.text.y = element_text(size=8,family = "Helvetica"),
        legend.text = element_text(size=8,family = "Helvetica"),
        legend.title = element_blank())
        #legend.justification = c(0.5, 0),
        #legend.position = "top")
pdf(file="../plots/numUTRsStaged.pdf",height=2.5,width=4.25,colormodel="rgb")
print(g2)
dev.off()

setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/novelUTRs/")
df <- read.table("novelUTRs.txt",header=TRUE,sep='\t')
# df$dataset <- factor(df$dataset,levels = c("Our genes","Waterston genes","Our isoforms","Waterston isoforms"))
g3 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(genes_with_novel_utrs,novel_utrs),labels=c("genes with novel UTRs"="genes with\nnovel UTRs"))+
  ylab("Number identified")+
  xlab("Stage")+
  theme(text = element_text(size = 8,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,family = "Helvetica"),
        axis.text.y = element_text(size=8,family = "Helvetica"),
        legend.text = element_text(size=8,family = "Helvetica"),
        legend.title = element_blank())
pdf(file="../plots/novelUTRsStaged.pdf",height=2.5,width=4.25,colormodel="rgb")
print(g3)
dev.off()
