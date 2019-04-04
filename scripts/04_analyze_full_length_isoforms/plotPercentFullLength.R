#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
# library(RColorBrewer)
library(extrafont)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!= 0){
fl_read_color = args[1]
}else{
  fl_read_color = '#389E34'
}

setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/plots/")
# COL <- (c(rgb(187,240,255,maxColorValue=255),rgb(119,226,255,maxColorValue=255),rgb(0,132,169,maxColorValue=255),rgb(0,99,127,maxColorValue=255),rgb(0,66,85,maxColorValue=255)))
# COL2 <- c("#FF050D","#FF050D7F","#8587FF","#8587FF7F","#00C52D50","#00C52D7F")
# COL_dRNA_cDNA <- c("#AAC78E","#AAC78E7F","#1F78B4","#1F78B47F")

# df1 <- data.frame(
#     sample= c(rep("L1_1",2),rep("L1_2",2),rep("L2_1",2),rep("L2_2",2),rep("L3_1",2),rep("L3_2",2),rep("L4_1",2),rep("L4_2",2),rep("YA_1",2),rep("YA_2",2),rep("ML_1",2)),
#     Metric=c(rep(c("Full Length","Short"),11)),
#     value = c(117929,71520,326284,181916,39394,24283,384922,202128,128576,73185,319404,171445,167140,101845,119149,70056,500411,311141,47994,29832,198883,136464) )

# p1 <- ggplot(df1,aes(fill=Metric,y=value,x=sample))+
#   geom_bar(stat="identity",position="fill")+
#   xlab("Stage")+
#   ylab("% Read Type")

# df1 <- data.frame(
#   sample = c("L1_1","L1_2","L2_1","L2_2","L3_1","L3_2","L4_1","L4_2","Young_Adult_1","Young_Adult_2","Gravid_Adult_1","Gravid_Adult_2","Male_1","Male_2"),
#   value = c(72.623954,
#             71.943968,
#             75.277491,
#             73.860038,
#             75.145717,
#             72.492229,
#             73.455176,
#             73.785229,
#             74.099695,
#             69.488068,
#             71.542334,
#             64.414973,
#             62.589102,
#             66.828008)
# )
# df1 <- data.frame(
#   sample = c("L1_1","L1_2","L2_1","L2_2","L3_1","L3_2","L4_1","L4_2","Young_Adult_1","Young_Adult_2","Gravid_Adult_1","Gravid_Adult_2","Male_1","Male_2"),
#   value = c(68.69511222,
#             67.42857143,
#             71.10965582,
#             69.6593988,
#             71.45393244,
#             68.71150464,
#             69.19783003,
#             70.10654905,
#             69.02986917,
#             66.72975069,
#             59.9798748,
#             63.4253738,
#             67.14481764,
#             61.51629623)
# )
# df1$sample <- factor(df1$sample,as.character(unique(df1$sample)))
# p1 <- ggplot(df1,aes(y=value,x=sample))+
#   geom_col(fill="#2579B2")+
#   ylim(c(0,100))+
#   ylab("% Full Length")+
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# 
# save_plot(filename = "full_length_percent.pdf",p1,base_aspect_ratio = 1.25)


df1 <- data.frame(
  sample = c("L1","L2","L3","L4","young\nadult ","mature\nadult ","male"),
  value = c(72.609,
            74.526,
            73.745,
            74.053,
            74.178,
            65.983,
            68.528),
  trimmed_value = c(72.6,
            74.5,
            73.7,
            74.1,
            74.2,
            66.0,
            68.5)
)
df1$sample <- factor(df1$sample,as.character(unique(df1$sample)))
p1 <- ggplot(df1,aes(y=value,x=sample))+
  geom_col(fill=fl_read_color,color="black",width=0.3)+
  geom_text(aes(label=trimmed_value),vjust=-1,size=3.5273)+
  ylab("% Full Length")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0),limits = c(0,100))+
  theme(text = element_text(size = 10,family="Helvetica"),
        axis.text.x = element_text(size=10,family="Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,family="Helvetica"))

pdf("full_length_percent_combined.pdf",height=2.5,width=4.167,colormodel="rgb")
p1
dev.off()
#save_plot(filename = "full_length_percent_combined.pdf",p1,base_aspect_ratio = 1.75)


#df2 <- data.frame(
#  sample= c("L1","L2","L3","L4","YA","ML","GA","All"),
#  Metric=c("Read Counts","Read Counts","Read Counts","Read Counts","Read Counts"),
#  value = c(697649,650727,692610,458190,889378,335347,727537,4451438)
#  )
# df2 <- data.frame(
#   sample= c(rep("L1",3),rep("L2",3),rep("L3",3),rep("L4",3),rep("YA",3),rep("ML",3),rep("GA",3),rep("All",3)),
#   Metric=c(rep(c("ExactMatch","WithinTolerance","NotWithinTolerance"),8)),
#   value = c(122432,2055,60980,
#             124349,2035,59083,
#             119894,2206,63367,
#             108244,2388,74835,
#             122216,1829,61422,
#             115977,2295,67195,
#             112745,1659,71063,
#             156909,1670,26888))
# positions <- c("L1","L2","L3","L4","YA","ML","GA","All")
#   
# p2 <- ggplot(df2,aes(fill=Metric,y=value,x=sample))+
#     geom_bar(stat="identity",position="fill")+
#     scale_x_discrete(limits=positions)+
#     xlab("Stage")+
#     ylab("% Splice Junctions")
#   
# save_plot(filename = "splice_quantification.pdf",p2,base_aspect_ratio = 2)

  
# df3 <- data.frame(
#   sample= c(rep("L1",4),rep("L2",4),rep("L3",4),rep("L4",4),rep("YA",4),rep("ML",4),rep("GA",4),rep("All",4)),
#   Metric=c(rep(c("ExactMatch","WithinTolerance","NotWithinTolerance","ReadCounts"),8)),
#   value = c(100*122432/185467,100*2055/185467,100*60980/185467,697649,
#             100*124349/185467,100*2035/185467,100*59083/185467,650727,
#             100*119894/185467,100*2206/185467,100*63367/185467,692610,
#             100*108244/185467,100*2388/185467,100*74835/185467,458190,
#             100*122216/185467,100*1829/185467,100*61422/185467,889378,
#             100*115977/185467,100*2295/185467,100*67195/185467,335347,
#             100*112745/185467,100*1659/185467,100*71063/185467,727537,
#             100*156909/185467,100*1670/185467,100*26888/185467,4451438)
# )
# df4 <- data.frame(
#   sample= c("L1","L2","L3","L4","YA","ML","GA","All"),
#   ExactMatches = c(100*122432/185467,100*124349/185467,100*119894/185467,100*108244/185467,100*122216/185467,100*115977/185467,100*112745/185467,100*156909/185467),
#   ReadCounts =c(697649,650727,692610,458190,889378,335347,727537,4451438)
# )
# 
# p3 <- ggplot(df4,aes(x=log(ReadCounts),y=ExactMatches))+
#   geom_point()+
#   ylab("Exact Matches %")
#p2 <- ggplot(df2,aes(factor(sample),value))+
#  geom_bar(stat="identity",position="dodge",alpha=0.8,fill=brewer.pal(5,"Set1")[3])+
#  scale_fill_brewer(palette="Greens")+
#  xlab("")+
# ylab("Number of Reads")+
#  scale_y_continuous(position="right")+
#  theme(text = element_text(size=20),axis.text = element_text(size=20))
#save_plot(filename = "readCounts.pdf",p2)

#df3 <- data.frame(
#  sample= c("L1","L2","L3","L4","YA"),
#  ReadCounts = c(711318,	665611,	708321,	469894,	912381),
#  Sensitivity = c(67.6,73.2,66.9,64.6,82.8) )

#p3 <-ggplot(df3,aes(x=ReadCounts,y=Sensitivity,color=sample))+
#  geom_point(size=4)
