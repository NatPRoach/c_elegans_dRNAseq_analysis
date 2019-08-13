#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(here)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!= 0){
fl_read_color = args[1]
}else{
  fl_read_color = '#389E34'
}

setwd(here())

# df1 <- data.frame(
#   sample = c("L1","L2","L3","L4","young\nadult ","mature\nadult ","male"),
#   value = c(72.609,
#             74.526,
#             73.745,
#             74.053,
#             74.178,
#             65.983,
#             68.528),
#   trimmed_value = c(72.6,
#             74.5,
#             73.7,
#             74.1,
#             74.2,
#             66.0,
#             68.5)
# )
# df1$sample <- factor(df1$sample,as.character(unique(df1$sample)))
# p1 <- ggplot(df1,aes(y=value,x=sample))+
#   geom_col(fill=fl_read_color,color="black",width=0.3)+
#   geom_text(aes(label=trimmed_value),vjust=-1,size=3.5273)+
#   ylab("% Full Length")+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0),limits = c(0,100))+
#   theme(text = element_text(size = 10,family="Helvetica"),
#         axis.text.x = element_text(size=10,family="Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size=10,family="Helvetica"))
# 
# pdf("figures/figure1/figure1C.pdf",height=2.5,width=4.167,colormodel="rgb")
# p1
# dev.off()


df1 <- data.frame(
  sample = c("L1","L2","L3","L4","young\nadult ","mature\nadult ","male","all"),
  value = c(59.97363612,
            62.33513332,
            62.63006464,
            61.50896704,
            56.8644868,
            57.70568005,
            51.88823248,
            58.75672198),
  trimmed_value = c(60.0,
                    62.3,
                    62.6,
                    61.5,
                    56.9,
                    57.7,
                    51.9,
                    58.8)
)
df1$sample <- factor(df1$sample,as.character(unique(df1$sample)))
p1 <- ggplot(df1,aes(y=value,x=sample))+
  geom_col(fill=fl_read_color,color="black",width=0.3)+
  geom_text(aes(label=trimmed_value),vjust=-1,size=3.5273)+
  ylab("% full-length")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0),limits = c(0,100))+
  theme(text = element_text(size = 10,family="Helvetica"),
        axis.text.x = element_text(size=10,family="Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10,family="Helvetica"))

# pdf("figures/figure1/stringent_figure1C.pdf",height=2.5,width=4.167,colormodel="rgb")
pdf("figures/figure1/figure1C.pdf",height=2.5,width=4.167,colormodel="rgb")
p1
dev.off()


# df1 <- data.frame(
#   sample = c("L1","L2","L3","L4","young\nadult ","mature\nadult ","male","all"),
#   value = c(72.192,
#             74.085,
#             73.255,
#             73.594,
#             73.889,
#             65.716,
#             67.638,
#             71.450),
#   trimmed_value = c(72.2,
#                     74.1,
#                     73.3,
#                     73.6,
#                     73.9,
#                     65.7,
#                     67.6,
#                     71.4)
# )
# df1$sample <- factor(df1$sample,as.character(unique(df1$sample)))
# p1 <- ggplot(df1,aes(y=value,x=sample))+
#   geom_col(fill=fl_read_color,color="black",width=0.3)+
#   geom_text(aes(label=trimmed_value),vjust=-1,size=3.5273)+
#   ylab("% Full Length")+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0),limits = c(0,100))+
#   theme(text = element_text(size = 10,family="Helvetica"),
#         axis.text.x = element_text(size=10,family="Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size=10,family="Helvetica"))
# 
# pdf("figures/figure1/sensitive_figure1C.pdf",height=2.5,width=4.167,colormodel="rgb")
# p1
# dev.off()