#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
# library(RColorBrewer)
library(here)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=0){
  color = args[1]
}else{
  color= "#C9B3D5"
}
setwd(here())

# df1 <- data.frame(
#   sample = c("AAUAAA","AltPAS","noPAS"),
#   value = c(100*5619/7732,100*4935/7559,100*283/1034)
# )
# 
# p1 <- ggplot(df1,aes(y=value,x=sample))+
#   geom_col(fill=color,color="black",width=0.3)+
#   ylim(c(0,100))+
#   ylab("% of Overlapping with Mangone")+
#   theme(text = element_text(size = 8,family="Helvetica"),
#         axis.text.x = element_text(size=8,family="Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size=8,family="Helvetica"))
# 
# pdf("figures/supplementals/sfigure4/sfigure4C.pdf",height=2.5,width=2.777,colormodel="rgb")
# p1
# dev.off()

# ### Sensitive
# df1 <- data.frame(
#   sample = c("AAUAAA","AltPAS","noPAS"),
#   value = c(100*5793/7898,100*5066/7767,100*281/1066)
# )
# 
# p1 <- ggplot(df1,aes(y=value,x=sample))+
#   geom_col(fill=color,color="black",width=0.3)+
#   ylim(c(0,100))+
#   ylab("% of Overlapping with Mangone")+
#   theme(text = element_text(size = 8,family="Helvetica"),
#         axis.text.x = element_text(size=8,family="Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size=8,family="Helvetica"))
# 
# pdf("figures/supplementals/sfigure4/sensitive_sfigure4C.pdf",height=2.5,width=2.777,colormodel="rgb")
# p1
# dev.off()


### Stringents
df1 <- data.frame(
  sample = c("AAUAAA","AltPAS","noPAS"),
  value = c(100*5496/7753,100*4855/7572,100*282/1017)
)

p1 <- ggplot(df1,aes(y=value,x=sample))+
  geom_col(fill=color,color="black",width=0.3)+
  ylim(c(0,100))+
  ylab("% of Overlapping with Mangone")+
  theme(text = element_text(size = 8,family="Helvetica"),
        axis.text.x = element_text(size=8,family="Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8,family="Helvetica"))

# pdf("figures/supplementals/sfigure4/stringent_sfigure4C.pdf",height=2.5,width=2.777,colormodel="rgb")
pdf("figures/supplementals/sfigure4/sfigure4C.pdf",height=2.5,width=2.777,colormodel="rgb")
p1
dev.off()