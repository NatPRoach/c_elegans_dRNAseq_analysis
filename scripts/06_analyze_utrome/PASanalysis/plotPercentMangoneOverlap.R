#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
# library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=0){
  color = args[1]
}else{
  color= "#C9B3D5"
}
setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/plots/")

df1 <- data.frame(
  sample = c("AAUAAA","AltPAS","noPAS"),
  value = c(100*5619/7732,100*4935/7559,100*283/1034)
)

p1 <- ggplot(df1,aes(y=value,x=sample))+
  geom_col(fill=color,color="black",width=0.3)+
  ylim(c(0,100))+
  ylab("% of Overlapping with Mangone")+
  theme(text = element_text(size = 8,family="Helvetica"),
        axis.text.x = element_text(size=8,family="Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8,family="Helvetica"))

pdf("percentMangoneOverlap.pdf",height=2.5,width=2.777,colormodel="rgb")
p1
dev.off()
