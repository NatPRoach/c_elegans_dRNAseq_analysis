#!/Library/Frameworks/R.framework/Resources/bin/Rscript

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
# library(RColorBrewer)
library(gdata)
library(extrafont)
library(here)
args = commandArgs(trailingOnly=TRUE)



if (length(args)!=0){
splice_isoform_color = args[1]
utr_color = args[2]
fl_isoform_color = args[3]
}else{
  splice_isoform_color= "#C9B3D5"
  utr_color = "#67A9CF"
  fl_isoform_color = "#02818A"
}

# #setwd("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countCombinatorialIsoforms/")
# 
# 
# 
# #+++++++++++++++++++++++++
# # Function to calculate the mean and the standard deviation
# # for each group
# #+++++++++++++++++++++++++
# # data : a data frame
# # varname : the name of a column containing the variable
# #to be summariezed
# # groupnames : vector of column names to be used as
# # grouping variables
# data_summary <- function(data, varname, groupnames){
#   require(plyr)
#   summary_func <- function(x, col){
#     c(mean = mean(x[[col]], na.rm=TRUE),
#       sd = sd(x[[col]], na.rm=TRUE))
#   }
#   data_sum<-ddply(data, groupnames, .fun=summary_func,
#                   varname)
#   data_sum <- rename(data_sum, c("mean" = varname))
#   return(data_sum)
# }
# 
# 
# 
# df <- read.table("combinatorial_isoform_count.txt",sep = "\t", header= TRUE)
# 
# df$dataset <- factor(df$dataset,levels = c("splice isoforms","utrs","full length isoforms"))
# g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
#   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_fill_manual(values=c(splice_isoform_color,utr_color,fl_isoform_color),labels = c("splice isoforms"="splice\nisoforms","utrs" = "3'UTRs","full length isoforms"="full length\nisoforms"))+
#   ylab("Number identified")+
#   xlab("Stage")+
#   theme(text = element_text(size = 8,family="Helvetica"),
#         axis.text.x = element_text(size=8,family="Helvetica"),
#         axis.text.y = element_text(size=8,family="Helvetica"))
#         # legend.title = element_blank(),
#         # legend.justification = c(0.5, 0),
#         # legend.position = "top")
# 
# 
# pdf(file="plots/numFullLengthIsoformsStaged.pdf",height=2.5,width=5,colormodel="rgb")
# print(g2)
# dev.off()
# 
# # df3 <- df[df$stage=="all",]
# # df3$dataset <- factor(df3$dataset,levels = c("splice isoforms","utrs","full length isoforms"))
# # g4 <- ggplot(df3,aes(x=dataset,y=counts))+
# #   geom_bar(position="dodge",stat = "identity",aes(fill =dataset),color="black",show.legend = FALSE)+
# #   #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
# #   scale_x_discrete(labels = c("utrs" = "3'UTRs"))+
# #   scale_fill_manual(values=c("#C9B3D5","#67A9CF","#02818A"),labels = c("utrs" = "3'UTRs"))+
# #   ylab("Number identified")+
# #   xlab("")+
# #   theme(axis.ticks.x=element_blank(),text = element_text(size = 20),axis.text.x = element_text(size=18),axis.text.y = element_text(size=20))
# # 
# # pdf(file="plots/numFullLengthIsoforms.pdf",height=6,width=8,colormodel="rgb")
# # print(g4)
# # dev.off()
# 
# df2 <- read.table("combinatorial_isoform_count_subsampled.txt",sep="\t",header=TRUE)
# # df2_summary <- data_summary(df2,varname = "counts",groupnames = c("stage","fraction_retained","dataset"))
# df2$stage <- factor(df2$stage,levels = c("L1","L2","L3","L4","young adult","mature adult","male","all"))
# 
# new <- data.frame(x = seq(1,3000000,10000))
# x<- df2$read_count[df2$stage=="L1"]
# y<- df2$counts[df2$stage=="L1"]
# # x <- x[51:100]
# # y <- y[51:100]
# #model <-nls(y~a/(1+exp(-b*(x-c))),start=c(a=15000,b=0.00001,c=166000))
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# L1 <- data.frame(read_count = new$x,counts=y_hat)
# l1_max <- max(x)
# l1_saturation <- coef(model)[[1]]
# 
# x<- df2$read_count[df2$stage=="L2"]
# y<- df2$counts[df2$stage=="L2"]
# # x <- x[51:100]
# # y <- y[51:100]
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# L2 <- data.frame(read_count = new$x,counts=y_hat)
# l2_max <- max(x)
# l2_saturation <- coef(model)[[1]]
# 
# x<- df2$read_count[df2$stage=="L3"]
# y<- df2$counts[df2$stage=="L3"]
# # x <- x[51:100]
# # y <- y[51:100]
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# L3 <- data.frame(read_count = new$x,counts=y_hat)
# l3_max <- max(x)
# l3_saturation <- coef(model)[[1]]
# 
# x<- df2$read_count[df2$stage=="L4"]
# y<- df2$counts[df2$stage=="L4"]
# # x <- x[51:100]
# # y <- y[51:100]
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# L4 <- data.frame(read_count = new$x,counts=y_hat)
# l4_max <- max(x)
# l4_saturation <- coef(model)[[1]]
# 
# x<- df2$read_count[df2$stage=="young adult"]
# y<- df2$counts[df2$stage=="young adult"]
# # x <- x[51:100]
# # y <- y[51:100]
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# young_adult <- data.frame(read_count = new$x,counts=y_hat)
# ya_max <- max(x)
# ya_saturation <- coef(model)[[1]]
# 
# 
# x<- df2$read_count[df2$stage=="mature adult"]
# y<- df2$counts[df2$stage=="mature adult"]
# # x <- x[51:100]
# # y <- y[51:100]
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# mature_adult <- data.frame(read_count = new$x,counts=y_hat)
# ga_max <- max(x)
# ga_saturation <- coef(model)[[1]]
# 
# x<- df2$read_count[df2$stage=="male"]
# y<- df2$counts[df2$stage=="male"]
# # x <- x[51:100]
# # y <- y[51:100]
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# male <- data.frame(read_count = new$x,counts=y_hat)
# ml_max <- max(x)
# ml_saturation <- coef(model)[[1]]
# 
# x<- df2$read_count[df2$stage=="all"]
# y<- df2$counts[df2$stage=="all"]
# # x <- x[51:100]
# # y <- y[51:100]
# model <- nls(y~a*x/(b+x),start=c(a=25000,b=100000))
# y_hat <- predict(model,new)
# all <- data.frame(read_count = new$x,counts=y_hat)
# 
# test <- combine(L1,L2,L3,L4,young_adult,mature_adult,male,all)
# test$stage <- factor(test$source,levels = c("L1","L2","L3","L4","young_adult","mature_adult","male","all"),
#                       labels= c("L1","L2","L3","L4","young adult","mature adult","male","all"))
# g3 <- ggplot(df2,aes(x=read_count,y=counts,color=stage))+
#   #geom_line(data = test,linetype=2)+
#   geom_point()+
#   #geom_errorbar(aes(ymin=counts-sd,ymax=counts+sd),width=0.05)+
#   scale_color_manual(values=c("#2579B2","#389E34","#E01E26","#FD7F23","#6A4098","#A7CEE2","#B3DE8E","#000000"))+
#   ylab("Number of full-length isoforms")+
#   xlab("Number of reads retained")+
#   theme(text = element_text(size = 10,family="Helvetica"),
#         axis.text.x = element_text(size=10,family="Helvetica"),
#         axis.text.y = element_text(size=10,family="Helvetica"))
# # g3 <- ggplot(df2_summary,aes(x=fraction_retained,y=counts,color=stage))+
# #   geom_line()+
# #   geom_point()+
# #   geom_errorbar(aes(ymin=counts-sd,ymax=counts+sd),width=0.05)+
# #   scale_color_manual(values=c("#2579B2","#389E34","#E01E26","#FD7F23","#6A4098","#A7CEE2","#B3DE8E","#F99B9B"))+
# #   ylab("Number of full-length isoforms")+
# #   xlab("Fraction of reads retained")
# 
# pdf(file="plots/numFullLengthIsoformsSubsampled.pdf",height=3,width=8,colormodel="rgb")
# print(g3)
# dev.off()
# 
# 
# L1 <- data.frame(max=l1_max,saturation=l1_saturation)
# L2 <- data.frame(max=l2_max,saturation=l2_saturation)
# L3 <- data.frame(max=l3_max,saturation=l3_saturation)
# L4 <- data.frame(max=l4_max,saturation=l4_saturation)
# YA <- data.frame(max=ya_max,saturation=ya_saturation)
# GA <- data.frame(max=ga_max,saturation=ga_saturation)
# ML <- data.frame(max=ml_max,saturation=ml_saturation)
# 
# All <- combine(L1,L2,L3,L4,YA,GA,ML)
# All$source <- factor(All$source,labels=c("L1"="L1",
#                                          "L2"="L2",
#                                          "L3"="L3",
#                                          "L4"="L4",
#                                          "YA"="young adult",
#                                          "GA" = "mature adult",
#                                          "ML"="male"))
# 
# g4 <- ggplot(All,aes(x=max,y=saturation))+
#   geom_point()+
#   geom_text(aes(label=source),nudge_y = 200,size=2.8219,family="Helvetica")+
#   xlab("Number of reads")+
#   ylab("Estimated saturation point")+
#   xlim(c(250000,550000))+
#   theme(text = element_text(size = 8,family="Helvetica"),
#         axis.text.x = element_text(size=8,family="Helvetica"),
#         axis.text.y = element_text(size=8,family="Helvetica"))
# pdf(file="plots/readsVsSaturation.pdf",height=3,width=5,colormodel="rgb")
# print(g4)
# dev.off()






