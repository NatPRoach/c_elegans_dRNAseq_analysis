#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(here)

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


setwd(here())

df <- read.table("results/scratch/PASanalysis/pasByStage.txt",sep="\t",header=TRUE)

l1 <- c(nrow(df[df$stage=="L1" & df$PAS_type=="AAUAAA",]),
        nrow(df[df$stage=="L1" & df$PAS_type=="altPAS",]),
        nrow(df[df$stage=="L1" & df$PAS_type=="noPAS",]))
l2 <- c(nrow(df[df$stage=="L2" & df$PAS_type=="AAUAAA",]),
        nrow(df[df$stage=="L2" & df$PAS_type=="altPAS",]),
        nrow(df[df$stage=="L2" & df$PAS_type=="noPAS",]))
l3 <- c(nrow(df[df$stage=="L3" & df$PAS_type=="AAUAAA",]),
        nrow(df[df$stage=="L3" & df$PAS_type=="altPAS",]),
        nrow(df[df$stage=="L3" & df$PAS_type=="noPAS",]))
l4 <- c(nrow(df[df$stage=="L4" & df$PAS_type=="AAUAAA",]),
        nrow(df[df$stage=="L4" & df$PAS_type=="altPAS",]),
        nrow(df[df$stage=="L4" & df$PAS_type=="noPAS",]))
ya <- c(nrow(df[df$stage=="young adult" & df$PAS_type=="AAUAAA",]),
        nrow(df[df$stage=="young adult" & df$PAS_type=="altPAS",]),
        nrow(df[df$stage=="young adult" & df$PAS_type=="noPAS",]))
ga <- c(nrow(df[df$stage=="mature adult" & df$PAS_type=="AAUAAA",]),
        nrow(df[df$stage=="mature adult" & df$PAS_type=="altPAS",]),
        nrow(df[df$stage=="mature adult" & df$PAS_type=="noPAS",]))
ml <- c(nrow(df[df$stage=="male" & df$PAS_type=="AAUAAA",]),
        nrow(df[df$stage=="male" & df$PAS_type=="altPAS",]),
        nrow(df[df$stage=="male" & df$PAS_type=="noPAS",]))
all <-c(nrow(df[df$PAS_type=="AAUAAA",]),
        nrow(df[df$PAS_type=="altPAS",]),
        nrow(df[df$PAS_type=="noPAS",]))

M <- as.table(rbind(l1,l2))
#fisher.test(M)
chisq.test(M)
M <- as.table(rbind(l2,l3))
#fisher.test(M)
chisq.test(M)
M <- as.table(rbind(l3,l4))
#fisher.test(M)
chisq.test(M)
M <- as.table(rbind(l4,ya))
#fisher.test(M)
chisq.test(M)
M <- as.table(rbind(ya,ga))
#fisher.test(M)
chisq.test(M)
M <- as.table(rbind(ya,ml))
#fisher.test(M)
chisq.test(M)
# 
# M <- as.table(rbind(l1,l2))
# fisher.test(M)
# chisq.test(M)
# M <- as.table(rbind(l1,l3))
# fisher.test(M)
# chisq.test(M)
# M <- as.table(rbind(l1,l4))
# fisher.test(M)
# chisq.test(M)
# M <- as.table(rbind(l1,ya))
# fisher.test(M)
# chisq.test(M)
# M <- as.table(rbind(l1,ga))
# fisher.test(M)
# chisq.test(M)
# M <- as.table(rbind(l1,ml))
# fisher.test(M)
# chisq.test(M)
# 
# M <- as.table(rbind(all,l1))
# chisq.test(M)
# M <- as.table(rbind(all,l2))
# chisq.test(M)
# M <- as.table(rbind(all,l3))
# chisq.test(M)
# M <- as.table(rbind(all,l4))
# chisq.test(M)
# M <- as.table(rbind(all,ya))
# chisq.test(M)
# M <- as.table(rbind(all,ga))
# chisq.test(M)
# M <- as.table(rbind(all,ml))
# chisq.test(M)

# chisq.test(l1,y=l2)
# chisq.test(l1,y=l3)
# chisq.test(l1,y=l4)
# chisq.test(l1,y=ya)
# chisq.test(l1,y=ga)
# chisq.test(l1,y=ml)
df2 <- df[df$stage!="Bartel" & df$stage != "Mangone",]
#g1 <-  ggplot(df,aes(reorder(stage,x_order)))+
g1 <-  ggplot(df2,aes(reorder(stage,x_order)))+
  geom_bar(aes(fill = PAS_type), position = position_fill(reverse = TRUE),color="black",width=0.5)+
  scale_fill_manual(values=c(canon_pas,alt_pas,no_pas),labels = c("AAUAAA" = "AAUAAA", "altPAS" = "alt PAS", "noPAS" = "no PAS"))+
  scale_y_continuous(labels = scales::percent_format())+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all","Mangone"= "Mangone","Bartel"="Jan"))+
  #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all","Mangone"= "Mangone","Bartel"="Jan"))+
  ylab("Percent PAS type")+
  xlab("dataset")+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
        #legend.key.size = unit(1.25, 'lines'),
        legend.position = "top",
        legend.title = element_blank(),
        #legend.justification = c(0.5, 0),
        legend.justification = c(0.5, 0))
pdf(file="figures/figure3/figure3F.pdf",height=2.5,width=4.1667,colormodel="rgb")
print(g1)
dev.off()

