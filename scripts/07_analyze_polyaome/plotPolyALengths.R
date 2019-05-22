#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(gdata)
library(here)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}


setwd(here())

tmp <- read.table("data/L1/bio1/tech1/analysis/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L1_1 <- tmp

tmp <- read.table("data/L1/bio1/tech2/analysis/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L1_2 <- tmp

tmp <- read.table("data/L2/bio1/tech1/analysis/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L2_1 <- tmp

tmp <- read.table("data/L2/bio1/tech2/analysis/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L2_2 <- tmp

tmp <- read.table("data/L3/bio1/tech1/analysis/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L3_1 <- tmp

tmp <- read.table("data/L3/bio1/tech2/analysis/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L3_2 <- tmp

tmp <- read.table("data/L4/bio1/tech1/analysis/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L4_1 <- tmp

tmp <- read.table("data/L4/bio1/tech2/analysis/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
L4_2 <- tmp

tmp <- read.table("data/male/bio1/tech1/analysis/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
ML_1 <- tmp

tmp <- read.table("data/male/bio1/tech2/analysis/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
ML_2 <- tmp

tmp <- read.table("data/young_adult/bio1/tech1/analysis/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
YA_1 <- tmp

tmp <- read.table("data/young_adult/bio1/tech2/analysis/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
YA_2 <- tmp


tmp <- read.table("data/adult/bio1/tech1/analysis/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
GA_1 <- tmp

tmp <- read.table("data/adult/bio1/tech2/analysis/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",header = TRUE,stringsAsFactors=FALSE)
tmp <- tmp[tmp$qc_tag == "PASS",]
GA_2 <- tmp

#merge by row
L1 <- rbind(L1_1,L1_2)
L2 <- rbind(L2_1,L2_2)
L3 <- rbind(L3_1,L3_2)
L4 <- rbind(L4_1,L4_2)
ML <- rbind(ML_1,ML_2)
YA <- rbind(YA_1,YA_2)
GA <- rbind(GA_1,GA_2)
all <- rbind(L1_1,L1_2,L2_1,L2_2,L3_1,L3_2,L4_1,L4_2,YA_1,YA_2,GA_1,GA_2,ML_1,ML_2)
All <- combine(L1,L2,L3,L4,YA,GA,ML,all)

larval <- combine(L1,L2,L3,L4)
adult <- combine(YA,GA,ML)

median(L1$polya_length)
median(L2$polya_length)
median(L3$polya_length)
median(L4$polya_length)
median(YA$polya_length)
median(GA$polya_length)
median(ML$polya_length)
median(larval$polya_length)
median(adult$polya_length)

p1 <- ggplot(All, aes(x = source, y = polya_length))+
  geom_hline(yintercept=25,color="gray60")+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=100,color="gray60")+
  geom_violin(aes(fill=source),show.legend = FALSE,draw_quantiles = c(0.5),color="black") +
  scale_x_discrete(labels=c("L1"="L1",
                            "L2"="L2",
                            "L3"="L3",
                            "L4"="L4",
                            "YA" = "yAd",
                            "GA" = "mAd",
                            "ML"= "male",
                            "all"="all"),expand=c(0.01,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,125))+
  scale_fill_manual(values=c(polya_color,polya_color,polya_color,polya_color,polya_color,polya_color,polya_color,polya_color))+
  xlab("Stage")+
  ylab("polyA tail lengths")+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=9,family="Helvetica",margin=margin(r=8,l=4)),
        legend.position = "top",
        legend.title = element_blank(),
        legend.justification = c(0.5, 0))

pdf(file="figures/figure4/figure4A.pdf",height=2.5,width=4.1667,colormodel="rgb")
print(p1)
dev.off()
