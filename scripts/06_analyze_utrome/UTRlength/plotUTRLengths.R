#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyr)
library(ggpubr)
library(here)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  our_utrs=args[1]
}else{
  our_utrs='#DFC27D'
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
# summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
#                       conf.interval=.95, .drop=TRUE) {
#   library(plyr)
#   
#   # New version of length which can handle NA's: if na.rm==T, don't count them
#   length2 <- function (x, na.rm=FALSE) {
#     if (na.rm) sum(!is.na(x))
#     else       length(x)
#   }
#   
#   # This does the summary. For each group's data frame, return a vector with
#   # N, mean, and sd
#   datac <- ddply(data, groupvars, .drop=.drop,
#                  .fun = function(xx, col) {
#                    c(N    = length2(xx[[col]], na.rm=na.rm),
#                      mean = mean   (xx[[col]], na.rm=na.rm),
#                      sd   = sd     (xx[[col]], na.rm=na.rm)
#                    )
#                  },
#                  measurevar
#   )
#   
#   # Rename the "mean" column    
#   datac <- rename(datac, c("mean" = measurevar))
#   
#   datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
#   
#   # Confidence interval multiplier for standard error
#   # Calculate t-statistic for confidence interval: 
#   # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
#   ciMult <- qt(conf.interval/2 + .5, datac$N-1)
#   datac$ci <- datac$se * ciMult
#   
#   return(datac)
# }
setwd(here())

df1 <- read.table("results/scratch/UTRlength/utr_lengths.txt",sep = "\t",header = TRUE )
# 
# median.quartile <- function(x){
#   out <- quantile(x, probs = c(0.25,0.5,0.75))
#   names(out) <- c("ymin","y","ymax")
#   return(out) 
# }


# df2 <- df1[df1$utr_length <= 500,]
# my_comparisons <- list( c("L1", "L2"), c("L2", "L3"), c("L3", "L4"), c("L4", "young adult"), c("young adult", "mature adult"), c("young adult", "male") )
# 
# g1 <- ggplot(df2,aes(x=reorder(stage,x_order),y=utr_length))+
#   geom_hline(yintercept=50,color="gray87")+
#   geom_hline(yintercept=150,color="gray87")+
#   geom_hline(yintercept=250,color="gray87")+
#   #stat_compare_means(comparisons = my_comparisons,label="p.signif")+
#   stat_compare_means(comparisons = my_comparisons,label="p.signif",label.y = c(520,540,520,540,520,560))+
#   # geom_violin(color="black",fill="#A7CEE2")+
#   geom_violin(draw_quantiles = c(0.50),color="black",fill="#A7CEE2")+
#   # stat_summary(fun.y=median.quartile,geom="point")+
#   ylim(c(0,575))+
#   ylab("3'UTR Length")+
#   scale_x_discrete(labels=c("L1","L2","L3","L4","young\nadult", "mature\nadult", "male","all"))+
#   xlab("")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())


#g1 <- ggplot(df1,aes(x=reorder(stage,x_order),y=utr_length))+
#g1 <- ggplot(df1,aes(x=reorder(stage,x_order),y=utr_length))

g1 <- ggplot(df1,aes(x=reorder(stage,x_order),y=utr_length))+
  geom_hline(yintercept=50,color="gray60")+
  geom_hline(yintercept=150,color="gray60")+
  geom_hline(yintercept=250,color="gray60")+
  # geom_violin(color="black",fill="#A7CEE2")+
  geom_violin(draw_quantiles = c(0.50),color="black",fill=our_utrs)+
  # stat_summary(fun.y=median.quartile,geom="point")+
  ylab("3'UTR Length")+
  #scale_x_discrete(labels=c("L1","L2","L3","L4","young\nadult", "mature\nadult", "male","all"))+
  scale_x_discrete(labels=c("L1","L2","L3","L4","yAd", "mAd", "male","all"))+
  xlab("")+
  expand_limits(y=0)+
  scale_y_continuous(limits=c(0,500),expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family = "Helvetica"),
        legend.title = element_blank())

pdf("figures/figure3/figure3E.pdf",height=2.5,width=4.25,colormodel="rgb")
g1
dev.off()
# 
# ks.test(df1$utr_length[df1$stage=="L1"],df1$utr_length[df1$stage=="L2"])
# ks.test(df1$utr_length[df1$stage=="L2"],df1$utr_length[df1$stage=="L3"])
# ks.test(df1$utr_length[df1$stage=="L3"],df1$utr_length[df1$stage=="L4"])
# ks.test(df1$utr_length[df1$stage=="L4"],df1$utr_length[df1$stage=="young adult"])
# ks.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="mature adult"])
# ks.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="male"])
# 
# wilcox.test(df1$utr_length[df1$stage=="L1"],df1$utr_length[df1$stage=="L2"])
# wilcox.test(df1$utr_length[df1$stage=="L2"],df1$utr_length[df1$stage=="L3"])
# wilcox.test(df1$utr_length[df1$stage=="L3"],df1$utr_length[df1$stage=="L4"])
# wilcox.test(df1$utr_length[df1$stage=="L4"],df1$utr_length[df1$stage=="young adult"])
# wilcox.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="mature adult"])
# wilcox.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="male"])


# g2 <- ggplot(df1,aes(x=utr_length,color=stage))+
#   #geom_hline(yintercept=100,color="gray87")+
#   #geom_hline(yintercept=150,color="gray87")+
#   #geom_hline(yintercept=200,color="gray87")+
#   stat_density(aes(group=stage),position="identity",geom="line")+
#   #ylim(c(0,600))+
#   xlim(0,1000)
#   xlab("3'UTR Length")+
#   #scale_x_discrete(labels=c("L1","L2","L3","L4","young\nadult", "mature\nadult", "male","all"))+
#   #xlab("")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())
#   
#tgc <- summarySE(df1,measurevar = "utr_length",groupvars=c("stage","x_order"))
# 
# g3 <- ggplot(tgc,aes(x=reorder(stage,x_order),y=utr_length))+
#   geom_bar(position = "dodge",stat="identity",fill="#2579B2",color="black",width = 0.3)+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_fill_manual(values=c("#2579B2"),labels=c(""))+
#   ylab("Average 3'UTR Length")+
#   xlab("Stage")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())
# pdf("../plots/averageUTRLengths.pdf",height=2.5,width=4.25,colormodel="rgb")
# g3
# dev.off()

# wilcox.test(df1$utr_length[df1$stage=="young adult"],df1$utr_length[df1$stage=="male"])
# wilcox.test(df1$utr_length[df1$stage=="L4"],df1$utr_length[df1$stage=="male"])
# 
# g1 <- ggplot(df1,aes(x=reorder(stage,x_order),y=utr_length))+
#   geom_violin(draw_quantiles = c(0.50))+
#   ylim(c(0,500))
# 
# g3 <- ggplot(df1)
# 
# df2 <- read.table("all_utr_lengths.txt",sep = "\t",header = TRUE )
# 
# g2 <- ggplot(df2,aes(x=reorder(stage,x_order),y=utr_length))+
#   geom_violin(draw_quantiles = c(0.5))+
#   ylim(c(0,500))
# 
# 
# tgc <- summarySE(df2,measurevar = "utr_length",groupvars=c("stage","x_order"))
# 
# g4 <- ggplot(tgc,aes(x=reorder(stage,x_order),y=utr_length))+
#   geom_bar(position = "dodge",stat="identity")
