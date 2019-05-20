#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(venn)
library(UpSetR)
library(eulerr)
library(here)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  mangone_utrs=args[1]
  jan_utrs=args[2]
  our_utrs=args[3]
}else{
  mangone_utrs='#B3DE8E'
  jan_utrs='#F99B9B'
  our_utrs='#A7CEE2'
}

setwd(here())

df <- read.table("results/overlaps/utr.overlap.matrix",sep="\t")
#colnames(df) <- c("Our UTRs","Mangone et al", "Jan et al","Wormbase")
colnames(df) <- c("This study","Mangone et al", "Jan et al")
e <- euler(df)

pdf(file="figures/figure3/figure3B.pdf",height=2.5,width=4.166,colormodel="rgb")
plot(e,fills = c(our_utrs,mangone_utrs,jan_utrs),labels=list(fontsize=8),quantities=list(fontsize=8))
dev.off()
# pdf(file="plots/UTRomeOverlaps.pdf",height=4.5,width=4.5,colormodel="rgb")
# venn(df,zcolor="style")
# dev.off()

# pdf(file="plots/UTRomeOverlapsUpSet.pdf",height=8,width=12,colormodel="rgb",onefile = F)
# #upset(df,order.by="freq",nsets=4,sets=c("Male","Gravid Adult", "Young Adult", "L4",  "L3", "L2","L1"),  point.size = 3.5, number.angles = 25, keep.order = TRUE,text.scale = c(3, 3, 1, 1.5, 3, 1.5))
# upset(df,order.by="freq",nsets=4, point.size = 3.5, number.angles = 25, keep.order = TRUE,text.scale = c(3, 3, 1, .8, 3, 1.5))
# #upset(df,order.by="freq",nsets=4, point.size = 3.5, number.angles = 25, keep.order = TRUE)#,text.scale = c(3, 3, 1, 1.5, 3, 1.5))
# 
# dev.off()