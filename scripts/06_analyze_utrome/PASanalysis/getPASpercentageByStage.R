#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  canon_pas=args[1]
  alt_pas=args[2]
  no_pas=args[3]
}else{
  canon_pas='#2579B2'
  alt_pas='#389E34'
  no_pas='#FD7F23'
}


setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/PASanalysis/")

df <- read.table("pasByStage.txt",sep="\t",header=TRUE)

g1 <-  ggplot(df,aes(reorder(stage,x_order)))+
  geom_bar(aes(fill = PAS_type), position = position_fill(reverse = TRUE),color="black",width=0.5)+
  scale_fill_manual(values=c(canon_pas,alt_pas,no_pas))+
  scale_y_continuous(labels = scales::percent_format())+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all","Mangone"= "Mangone","Bartel"="Jan"))+
  ylab("Percent PAS type")+
  xlab("dataset")+
  theme(text = element_text(size = 8,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,family = "Helvetica"),
        axis.text.y = element_text(size=8,family = "Helvetica"),
        legend.text = element_text(size=8,family = "Helvetica"),
        legend.title = element_blank(),
        legend.justification = c(0.5, 0),
        legend.position = "top")
pdf(file="../plots/pasStaged.pdf",height=2.5,width=4.5,colormodel="rgb")
print(g1)
dev.off()

