#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
#library(ggridges)
#library(lubridate)
#library(colorspace)
#library("dviz.supp")
library(gdata)
library(here)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  polya_color=args[1]
}else{
  polya_color='#A7CEE2'
}


setwd(here())

# theme_minimal_grid <- function(font_size = 14, font_family = "", line_size = .5,
#                                rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14,
#                                colour = "grey90") {
#   # Starts with theme_cowplot and then modifies some parts
#   theme_cowplot(font_size = font_size, font_family = font_family, line_size = line_size,
#                 rel_small = rel_small, rel_tiny = rel_tiny, rel_large = rel_large) %+replace%
#     theme(
#       # make grid lines
#       panel.grid        = element_line(colour = colour,
#                                        size = line_size),
#       panel.grid.minor  = element_blank(),
#       
#       # adjust axis tickmarks
#       axis.ticks        = element_line(colour = colour, size = line_size),
#       
#       # no x or y axis lines
#       axis.line.x       = element_blank(),
#       axis.line.y       = element_blank(),
#       
#       # no filled background for facted plots
#       strip.background = element_blank(),
#       
#       complete = TRUE
#     )
# }
# 
# 
# 
# theme_dviz_grid <- function(font_size = 14, font_family = dviz_font_family, line_size = .5,
#                             rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14,
#                             colour = "grey90") {
#   half_line <- font_size / 2
#   
#   theme_minimal_grid(font_size = font_size, font_family = font_family, line_size = line_size,
#                               rel_small = rel_small, rel_tiny = rel_tiny, rel_large = rel_large,
#                               colour = colour)  %+replace%
#     theme(
#       plot.margin = margin(half_line/2, 1.5, half_line/2, 1.5),
#       complete = TRUE
#     )
# }

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
# All <- combine(L1_1,L1_2,L2_1,L2_2,L3_1,L3_2,L4_1,L4_2,YA_1,YA_2,GA_1,GA_2,ML_1,ML_2)

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
# test <- ks.test(larval$polya_length,adult$polya_length)
# test <- wilcox.test(larval$polya_length,adult$polya_length)

# bandwidth <- 3.4
# polya_base <- ggplot(data=All, aes(x = polya_length,y = source,fill = ..x..))+
#   geom_density_ridges_gradient(
#     scale = 3, rel_min_height = 0.01, bandwidth = bandwidth,
#     color="black", size = 0.25,quantile_lines=TRUE,quantiles= 2
#   )+
#   scale_x_continuous(
#     name = "Estimated polyA length",
#     expand = c(0, 0), breaks = c(0, 25, 50, 75,100,125,150,175), labels = c("0","25","50","75","100","125","150","175"),limits=c(0,175)
#   ) +
#   theme_set(theme_cowplot(font_size=20))+
#   scale_y_discrete(labels=c("L1","L2","L3","L4","young adult", "mature adult", "male"),name = NULL, expand = c(0, .2, 0, 2.6)) +
#   scale_fill_continuous_sequential(palette = "Heat",p2=2, l1 = 20, l2 = 100, c2 = 0) +
#   guides(fill = "none") +
#   #theme_dviz_grid() +
#   theme(
#     axis.text.y = element_text(vjust = 0),
#     plot.margin = margin(3, 7, 3, 1.5)
#   )
# 
# pdf("plots/polyAdt_ridges.pdf")
# polya_base
# dev.off()
# 
# wilcox.test(L1$polya_length,L2$polya_length)
# wilcox.test(L2$polya_length,L3$polya_length)
# wilcox.test(L3$polya_length,L4$polya_length)
# wilcox.test(L4$polya_length,YA$polya_length)
# wilcox.test(YA$polya_length,GA$polya_length)
# wilcox.test(YA$polya_length,ML$polya_length)

# kruskal.test(list(L1$polya_length,L2$polya_length))$p.value
# kruskal.test(list(L2$polya_length,L3$polya_length))$p.value
# kruskal.test(list(L3$polya_length,L4$polya_length))$p.value
# kruskal.test(list(L4$polya_length,YA$polya_length))$p.value
# kruskal.test(list(YA$polya_length,GA$polya_length))$p.value
# kruskal.test(list(YA$polya_length,ML$polya_length))$p.value



# ks.test(L1$polya_length,L2$polya_length)
# ks.test(L2$polya_length,L3$polya_length)
# ks.test(L3$polya_length,L4$polya_length)
# ks.test(L4$polya_length,YA$polya_length)
# ks.test(YA$polya_length,GA$polya_length)
# ks.test(YA$polya_length,ML$polya_length)


# 
# L1_no_outliers <- L1[!L1$polya_length %in% boxplot.stats(L1$polya_length)$out,]
# L2_no_outliers <- L2[!L2$polya_length %in% boxplot.stats(L2$polya_length)$out,]
# L3_no_outliers <- L3[!L3$polya_length %in% boxplot.stats(L3$polya_length)$out,]
# L4_no_outliers <- L4[!L4$polya_length %in% boxplot.stats(L4$polya_length)$out,]
# YA_no_outliers <- YA[!YA$polya_length %in% boxplot.stats(YA$polya_length)$out,]
# GA_no_outliers <- GA[!GA$polya_length %in% boxplot.stats(GA$polya_length)$out,]
# ML_no_outliers <- ML[!ML$polya_length %in% boxplot.stats(ML$polya_length)$out,]
# All_no_outliers <- combine(L1_no_outliers,L2_no_outliers,L3_no_outliers,L4_no_outliers,YA_no_outliers,GA_no_outliers,ML_no_outliers)
# 
# ks.test(L1_no_outliers$polya_length,L2_no_outliers$polya_length,exact = T)
# ks.test(L2_no_outliers$polya_length,L3_no_outliers$polya_length)
# ks.test(L3_no_outliers$polya_length,L4_no_outliers$polya_length)
# ks.test(L4_no_outliers$polya_length,YA_no_outliers$polya_length)
# ks.test(YA_no_outliers$polya_length,GA_no_outliers$polya_length)
# ks.test(YA_no_outliers$polya_length,ML_no_outliers$polya_length)
# 
# 
# wilcox.test(L1_no_outliers$polya_length,L2_no_outliers$polya_length)
# wilcox.test(L2_no_outliers$polya_length,L3_no_outliers$polya_length)
# wilcox.test(L3_no_outliers$polya_length,L4_no_outliers$polya_length)
# wilcox.test(L4_no_outliers$polya_length,YA_no_outliers$polya_length)
# wilcox.test(YA_no_outliers$polya_length,GA_no_outliers$polya_length)
# wilcox.test(YA_no_outliers$polya_length,ML_no_outliers$polya_length)
# 
# L1_subsampled = L1_no_outliers[sample(nrow(L1_no_outliers), 10000), ]
# L2_subsampled = L2_no_outliers[sample(nrow(L2_no_outliers), 10000), ]
# L3_subsampled = L3_no_outliers[sample(nrow(L3_no_outliers), 10000), ]
# L4_subsampled = L4_no_outliers[sample(nrow(L4_no_outliers), 10000), ]
# YA_subsampled = YA_no_outliers[sample(nrow(YA_no_outliers), 10000), ]
# GA_subsampled = GA_no_outliers[sample(nrow(GA_no_outliers), 10000), ]
# ML_subsampled = ML_no_outliers[sample(nrow(ML_no_outliers), 10000), ]
# All_subsampled <- combine(L1_subsampled,L2_subsampled,L3_subsampled,L4_subsampled,YA_subsampled,GA_subsampled,ML_subsampled)
# 
# #L1 -> L2 shift
# wilcox.test(L1_subsampled$polya_length,L2_subsampled$polya_length)
# wilcox.test(L1_subsampled$polya_length,L3_subsampled$polya_length)
# wilcox.test(L1_subsampled$polya_length,L4_subsampled$polya_length)
# 
# #L2 - L4
# wilcox.test(L2_subsampled$polya_length,L3_subsampled$polya_length)
# wilcox.test(L3_subsampled$polya_length,L4_subsampled$polya_length)
# wilcox.test(L2_subsampled$polya_length,L4_subsampled$polya_length)
# 
# 
# 
# #L4 -> Adult stages shifts
# wilcox.test(L4_subsampled$polya_length,YA_subsampled$polya_length)
# 
# #Adult stages
# wilcox.test(YA_subsampled$polya_length,GA_subsampled$polya_length)
# wilcox.test(YA_subsampled$polya_length,ML_subsampled$polya_length)
# wilcox.test(ML_subsampled$polya_length,GA_subsampled$polya_length)
# 
p1 <- ggplot(All, aes(x = source, y = polya_length))+
  geom_hline(yintercept=25,color="gray60")+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=100,color="gray60")+
  # geom_hline(yintercept=40,color="gray60")+
  # geom_hline(yintercept=80,color="gray60")+
  # geom_hline(yintercept=120,color="gray60")+
  #theme(axis.text.x = element_text(angle=45,hjust=1))+
  geom_violin(aes(fill=source),show.legend = FALSE,draw_quantiles = c(0.5),color="black") +
  #geom_boxplot(aes(fill=source),show.legend = FALSE, color="black")+
  scale_x_discrete(labels=c("L1"="L1",
                            "L2"="L2",
                            "L3"="L3",
                            "L4"="L4",
                            "YA" = "yAd",
                            "GA" = "mAd",
                            "ML"= "male",
                            "all"="all"),expand=c(0.01,0))+
  # scale_x_discrete(labels=c("L1"="L1",
  #                           "L2"="L2",
  #                           "L3"="L3",
  #                           "L4"="L4",
  #                           "YA" = "young\nadult",
  #                           "GA" = "mature\nadult",
  #                           "ML"= "male",
  #                           "all"="all"),expand=c(0.01,0))+
  #scale_y_continuous(expand=c(0,0),limits=c(0,125))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(0,125))+
  #scale_color_manual(values=c("#a7cee2","#b3de8e","#b3de8e","#b3de8e","#f99b9b","#f99b9b","#f99b9b"))+
  scale_fill_manual(values=c(polya_color,polya_color,polya_color,polya_color,polya_color,polya_color,polya_color,polya_color))+
  xlab("Stage")+
  ylab("polyA tail lengths")+
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
#legend.justification = c(0.5, 0),
#legend.position = "top")

pdf(file="figures/figure4/figure4A.pdf",height=2.5,width=4.1667,colormodel="rgb")
print(p1)
dev.off()

# p1 <- ggplot(All_no_outliers, aes(x = source, y = polya_length))+
#   geom_hline(yintercept=40,color="gray87")+
#   geom_hline(yintercept=80,color="gray87")+
#   geom_hline(yintercept=120,color="gray87")+
#   theme_set(theme_cowplot(font_size=20))+
#   theme(axis.text.x = element_text(angle=45,hjust=1))+
#   geom_violin(aes(fill=source),show.legend = FALSE,draw_quantiles = c(0.5)) +
#   scale_x_discrete(labels=c("L1","L2","L3","L4","young adult", "mature adult", "male"),expand=c(0.01,0))+
#   scale_y_continuous(expand=c(0,0),limits=c(0,125))+
#   scale_color_manual(values=c("#a7cee2","#b3de8e","#b3de8e","#b3de8e","#f99b9b","#f99b9b","#f99b9b"))+
#   scale_fill_manual(values=c("#a7cee2","#b3de8e","#b3de8e","#b3de8e","#f99b9b","#f99b9b","#f99b9b"))+
#   xlab("Stage")+
#   ylab("polyA tail lengths")
# 
# pdf(file="polyAdt_no_outliers_colored.pdf",height=4.5,width=6,colormodel="rgb")
# print(p1)
# dev.off()
# 
# 
# p1 <- ggplot(All_no_outliers, aes(x = source, y = polya_length))+
#   geom_hline(yintercept=40,color="gray87")+
#   geom_hline(yintercept=80,color="gray87")+
#   geom_hline(yintercept=120,color="gray87")+
#   theme_set(theme_cowplot(font_size=20))+
#   theme(axis.text.x = element_text(angle=45,hjust=1))+
#   geom_violin(aes(fill=source),show.legend = FALSE,draw_quantiles = c(0.5)) +
#   scale_x_discrete(labels=c("L1","L2","L3","L4","young adult", "mature adult", "male"),expand=c(0.01,0))+
#   scale_y_continuous(expand=c(0,0),limits=c(0,125))+
#   scale_color_manual(values=c("#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2"))+
#   scale_fill_manual(values=c("#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2"))+
#   xlab("Stage")+
#   ylab("polyA tail lengths")
# 
# pdf(file="polyAdt_no_outliers.pdf",height=4.5,width=6,colormodel="rgb")
# print(p1)
# dev.off()
# 
# 
# # medians = data.frame(medians = c(median(L1$polya_length),
# #                                  median(L2$polya_length),
# #                                  median(L3$polya_length),
# #                                  median(L4$polya_length),
# #                                  median(YA$polya_length),
# #                                  median(GA$polya_length),
# #                                  median(ML$polya_length)),
# #                     source = c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),
# #                     order = c(0,1,2,3,4,5,6))
# # p1 <- ggplot(medians,aes(x=reorder(source,order),y=medians,fill="test"))+
# #   geom_bar(stat="identity",show.legend = FALSE)+
# #   scale_fill_manual(values=c("#2579B2"))+
# #   ylab("Median Poly(A) length")+
# #   theme(text = element_text(size=20),axis.text.x = element_text(size= 20, angle=30,hjust= 1),legend.title=element_blank(),axis.title.x=element_blank())
# 
#   
# # t.test(L1_no_outliers$polya_length,L2_no_outliers$polya_length)
# # t.test(L1_no_outliers$polya_length,L3_no_outliers$polya_length)
# # t.test(L1_no_outliers$polya_length,L4_no_outliers$polya_length)
# # t.test(L1_no_outliers$polya_length,YA_no_outliers$polya_length)
# # t.test(L1_no_outliers$polya_length,GA_no_outliers$polya_length)
# # t.test(L1_no_outliers$polya_length,ML_no_outliers$polya_length)
# # 
# # t.test(L3_no_outliers$polya_length,L4_no_outliers$polya_length)
# # 
# # t.test(YA_no_outliers$polya_length,GA_no_outliers$polya_length)
# # 
# # 
# # t.test(L1$polya_length,L2$polya_length)
# # t.test(L1$polya_length,L3$polya_length)
# # t.test(L1$polya_length,L4$polya_length)
# # t.test(L1$polya_length,YA$polya_length)
# # t.test(L1$polya_length,GA$polya_length)
# # t.test(L1$polya_length,ML$polya_length)
# # 
# # t.test(L2$polya_length,L3$polya_length)
# # t.test(L2$polya_length,L4$polya_length)
# # t.test(L2$polya_length,YA$polya_length)
# # t.test(L2$polya_length,GA$polya_length)
# # t.test(L2$polya_length,ML$polya_length)
# # 
# # t.test(L3$polya_length,L4$polya_length)
# # 
# # t.test(YA$polya_length,GA$polya_length)
# # t.test(YA$polya_length,ML$polya_length)
# 
# p1 <- ggplot(All, aes(x = source, y = polya_length))+
#   theme_set(theme_cowplot(font_size=20))+
#   theme(axis.text.x = element_text(angle=45,hjust=1))+
#   geom_violin(aes(fill=source),show.legend = FALSE,draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_x_discrete(labels=c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),expand=c(0.01,0))+
#   scale_y_continuous(expand=c(0,0),limits=c(0,125))+
#   scale_fill_manual(values=c("#1abcc2","#1abcc2","#1abcc2","#1abcc2","#1abcc2","#1abcc2","#1abcc2"))+
#   xlab("Stage")+
#   ylab("polyA tail lengths")
# 
# # p1 <- ggplot(All, aes(x = source, y = polya_length))+
# #   theme_set(theme_cowplot(font_size=20))+
# #   theme(axis.text.x = element_text(angle=45,hjust=1))+
# #   geom_violin(aes(fill=source),show.legend = FALSE) +
# #   geom_boxplot(width=0.2,aes(fill=source,color=source),show.legend = FALSE) +
# #   scale_x_discrete(labels=c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),expand=c(0.01,0))+
# #   scale_y_continuous(expand=c(0,0),limits=c(0,125))+
# #   scale_color_manual(values=c("#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2"))+
# #   scale_fill_manual(values=c("#1f6a98","#1f6a98","#1f6a98","#1f6a98","#1f6a98","#1f6a98","#1f6a98"))+
# #   xlab("Stage")+
# #   ylab("polyA tail lengths")
# 
# p1 <- ggplot(All, aes(x = source, y = polya_length))+
#   geom_hline(yintercept=40,color="gray87")+
#   geom_hline(yintercept=80,color="gray87")+
#   geom_hline(yintercept=120,color="gray87")+
#   theme_set(theme_cowplot(font_size=20))+
#   theme(axis.text.x = element_text(angle=45,hjust=1))+
#   geom_violin(aes(fill=source),show.legend = FALSE,draw_quantiles = c(0.5)) +
#   scale_x_discrete(labels=c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),expand=c(0.01,0))+
#   scale_y_continuous(expand=c(0,0),limits=c(0,125))+
#   scale_color_manual(values=c("#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2"))+
#   scale_fill_manual(values=c("#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2"))+
#   xlab("Stage")+
#   ylab("polyA tail lengths")
# 
# p2 <- ggplot(All, aes(x = source, y = polya_length))+
#   geom_hline(yintercept=40,color="gray87")+
#   geom_hline(yintercept=80,color="gray87")+
#   geom_hline(yintercept=120,color="gray87")+
#   theme_set(theme_cowplot(font_size=20))+
#   theme(axis.text.x = element_text(angle=45,hjust=1))+
#   geom_boxplot(aes(fill=source),show.legend = FALSE) +
#   scale_x_discrete(labels=c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),expand=c(0.01,0))+
#   scale_y_continuous(expand=c(0,0),limits=c(0,125))+
#   scale_fill_manual(values=c("#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2","#a7cee2"))+
#   xlab("Stage")+
#   ylab("polyA tail lengths")
# 
# plot_grid(
#   p1,
#   p2,
#   ncol = 1,
#   align = "hv",
#   rel_heights = c(6,9))
# 
# # p1 <- ggplot(All_no_outliers, aes(x = source, y = polya_length))+
# #   theme_set(theme_cowplot(font_size=20))+
# #   theme(axis.text.x = element_text(angle=45,hjust=1))+
# #   geom_boxplot(aes(fill=source),show.legend = FALSE) +
# #   scale_x_discrete(labels=c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),expand=c(0.01,0))+
# #   scale_y_continuous(expand=c(0,0),limits=c(0,200))+
# #   xlab("Stage")+
# #   ylab("polyA tail lengths")
# 
# 
# 
# 
# # p1 <- ggplot(All_no_outliers, aes(x = source, y = polya_length))+
# #   theme_set(theme_cowplot(font_size=20))+
# #   theme(axis.text.x = element_text(angle=45,hjust=1))+
# #   geom_violin(aes(fill=source),show.legend = FALSE) +
# #   scale_x_discrete(labels=c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),expand=c(0.01,0))+
# #   scale_y_continuous(expand=c(0,0),limits=c(0,200))+
# #   xlab("Stage")+
# #   ylab("polyA tail lengths")
# # p1 <- ggplot(All, aes(x = source, y = polya_length))+
# #   theme_set(theme_cowplot(font_size=20))+
# #   theme(axis.text.x = element_text(angle=45,hjust=1))+
# #   geom_violin(aes(fill=source),show.legend = FALSE) +
# #   scale_x_discrete(labels=c("L1","L2","L3","L4","Young Adult", "Gravid Adult", "Male"),expand=c(0.01,0))+
# #   scale_y_continuous(expand=c(0,0),limits=c(0,200))+
# #   xlab("Stage")+
# #   ylab("polyA tail lengths")
# 
# # p1 <- ggplot(All, aes(x = polya_length, y = source,fill=source))+
# #   theme_set(theme_cowplot(font_size=18))+
# #   guides(fill=FALSE)+
# #   geom_density_ridges(scale = 0.9) +
# #   scale_y_discrete(expand=c(0.01,0))+
# #   scale_x_continuous(expand=c(0,0),limits=c(0,200))+
# #   ylab("Stage")
# 
# # p1 <- ggplot(All, aes(x = polya_length, y = source))+
# #   theme_set(theme_cowplot(font_size=40))+
# #   guides(fill=FALSE)+
# #   geom_density_ridges(scale = 0.9) +
# #   scale_y_discrete(expand=c(0.01,0))+
# #   scale_x_continuous(expand=c(0,0),limits=c(0,200))+
# #   ylab("Stage")
# pdf(file="polyAdt_large_text2.pdf",height=4.5,width=6,colormodel="rgb")
# print(p1)
# dev.off()
# pdf(file="polyAdt_large_text3.pdf",height=4.5,width=6,colormodel="rgb")
# print(p2)
# dev.off()
# 
# #p2 <- ggplot(All,aes(x=polya_length))+
# #  geom_histogram(binwidth = 1, aes(y=..density..))+
# #  geom_density(alpha=0.2,fill="#FF6666")+
# #  scale_x_continuous(limits=c(0,200))
# p2 <- ggplot(All,aes(x=polya_length,fill="#a7cee2"))+
#   geom_histogram(binwidth = 1, aes(y=..density..),show.legend = FALSE)+
#   scale_fill_manual(values=c("#a7cee2"))+
#   scale_x_continuous(limits=c(0,200))
# pdf(file="polyAAll.pdf",height=4.5,width=6,colormodel="rgb")
# print(p2)
# dev.off()
# 
# p3 <- ggplot(L4,aes(x=polya_length,fill="#a7cee2"))+
#   geom_histogram(binwidth = 1, aes(y=..density..),show.legend = FALSE)+
#   scale_fill_manual(values=c("#a7cee2"))+
#   scale_x_continuous(limits=c(0,200))
# 
# pdf(file="polyAL4.pdf",height=4.5,width=6,colormodel="rgb")
# print(p3)
# dev.off()