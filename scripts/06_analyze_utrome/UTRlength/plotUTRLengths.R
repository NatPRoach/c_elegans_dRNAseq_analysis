#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(plyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=0){
  our_utrs=args[1]
}else{
  our_utrs='#A7CEE2'
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/UTRlength/")

df1 <- read.table("utr_lengths.txt",sep = "\t",header = TRUE )

# g1 <- ggplot(df1,aes(x=reorder(stage,x_order),y=utr_length))+
#   geom_violin(draw_quantiles = c(0.25,0.50,0.75))+
#   ylim(c(0,500))
tgc <- summarySE(df1,measurevar = "utr_length",groupvars=c("stage","x_order"))

g3 <- ggplot(tgc,aes(x=reorder(stage,x_order),y=utr_length))+
  geom_bar(position = "dodge",stat="identity",fill="#2579B2",color="black",width = 0.3)+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_fill_manual(values=c("#2579B2"),labels=c(""))+
  ylab("Average 3'UTR Length")+
  xlab("Stage")+
  theme(text = element_text(size = 8,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,family = "Helvetica"),
        axis.text.y = element_text(size=8,family = "Helvetica"),
        legend.text = element_text(size=8,family = "Helvetica"),
        legend.title = element_blank())
pdf("../plots/averageUTRLengths.pdf",height=2.5,width=4.25,colormodel="rgb")
g3
dev.off()

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
