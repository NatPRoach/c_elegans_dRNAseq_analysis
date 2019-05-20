#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
#library(grid)
#library(gridExtra)
#library(extrafont)
library(here)


args = commandArgs(trailingOnly=TRUE)
if (length(args)!=0){
all_read_color = args[1]
ful_read_color = args[2]
nfl_read_color = args[3]
gene_model_color = args[4]
}else{
  all_read_color = "#E01E26"
  ful_read_color = "#389E34"
  nfl_read_color = "#2579B2"
  gene_model_color = "#A7CEE2"
}

#setwd("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/metagene/data/")

#---L1-1--

all_data <-   read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.all_meta.txt"),header=FALSE)
full_data <-  read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.full_meta.txt"),header=FALSE)
short_data <- read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.short_meta.txt"),header=FALSE)

df <- data.frame(x = 1:200000,
                 all_reads = all_data,
                 full_reads = full_data,
                 short_reads = short_data)
colnames(df) <- c("x","all_reads","full_reads","short_reads")

# p1 <- ggplot(df) +
#   geom_line(aes(x=x, y=all_reads,color = "all reads"),size = 2)+
#   geom_line(aes(x=x, y=full_reads,color = "full-length\nreads"),size = 2)+
#   geom_line(aes(x=x, y=short_reads,color = "non full-length\nreads"),size = 2)+
#   geom_vline(xintercept=50000)+
#   geom_vline(xintercept=150000)+
#   # geom_text(data=data.frame(x=1:200),x=25000,  y=-40000,label = "5'UTR",family="Helvetica",size=8)+
#   # geom_text(data=data.frame(x=1:200),x=100000, y=-40000,label = "CDS",family="Helvetica",size=8)+
#   # geom_text(data=data.frame(x=1:200),x=175000, y=-40000,label = "3'UTR",family="Helvetica",size=8)+
#   # coord_cartesian( clip = 'off') +   # This keeps the labels from disappearing
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x= element_blank(),
#         axis.text.y = element_text(size=8,family="Helvetica"),
#         axis.title.y= element_text(size=8,family="Helvetica"),
#         legend.title = element_blank(),
#         legend.text = element_text(size=8,family="Helvetica"),
#         legend.key.size = unit(1.5, 'lines'))+
#   # xlab("")+
#   ylab("Relative Coverage")+
#   #scale_x_continuous(breaks= c(50000,150000),labels = c("CDS Start", "CDS Stop"))
#   scale_x_continuous(breaks= c(50000,150000),labels = c("", ""))+
#   scale_color_manual(values=c(all_read_color,ful_read_color,nfl_read_color))


p4 <- ggplot(df) +
  geom_vline(xintercept=50000,size=0.5)+
  geom_vline(xintercept=150000,size=0.5)+
  geom_line(aes(x=x, y=all_reads,color = "all reads"),size = 1)+
  geom_line(aes(x=x, y=full_reads,color = "full-length\nreads"),size = 1)+
  geom_line(aes(x=x, y=short_reads,color = "non full-length\nreads"),size = 1)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.y = element_text(size=10,family="Helvetica"),
        axis.title.y= element_text(size=10,family="Helvetica"),
        legend.title = element_blank(),
        legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
        legend.key.size = unit(1.25, 'lines'),
        legend.position = "top",
        #legend.justification = c(0.5, 0),
        legend.justification = c(0.5, 0))+
        #legend.spacing.x = unit(0.5, 'cm'))+
  geom_rect(data=data.frame(),aes(xmin=1, xmax=200000),ymin=-2.5e-7 - 3.5e-7- 1e-7,ymax=2.5e-7 - 3.5e-7- 1e-7,fill=gene_model_color,size =0.5,color="black")+
  geom_rect(data=data.frame(),aes(xmin=50000,xmax=150000),ymin=-3.5e-7 - 3.5e-7 - 2e-7,  ymax=3.5e-7 - 3.5e-7,  fill=gene_model_color,size =0.5,color="black")+
  # geom_segment(data=data.frame(),x=75000,xend=75000,y=0,yend=4,color="black",size = 1)+
  # geom_segment(data=data.frame(),x=110000,xend=110000,y=0,yend=4,color="black",size = 1)+
  geom_text(data=data.frame(x=1),x=25000, y=-7e-7 - 3.5e-7- 3e-7,label = "5'UTR",family="Helvetica",size=3.52734)+
  geom_text(data=data.frame(x=1),x=100000, y=-7e-7 - 3.5e-7- 3e-7,label = "CDS",family="Helvetica",size=3.52734)+
  geom_text(data=data.frame(x=1),x=175000, y=-7e-7 - 3.5e-7- 3e-7,label = "3'UTR",family="Helvetica",size=3.52734)+
  #geom_segment(x=50000,xend=50000,y=0,yend=4,color="black",size = 1)+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(-12e-7- 2e-7,8.25e-6),clip="off")+
  #geom_segment(x=50000,xend=50000,y=2,yend=6,color="black",size = 2,lineend="square")+
  #geom_segment(x=50000,xend=75000,y=6,yend=6, color="black",arrow=arrow(),size=2)+
  #geom_segment(x=50000,xend=51000,y=6,yend=6, color="black",size=2,lineend="square")+
  #blank_axes_and_thin_margin
  ylab("Relative Coverage")+
  scale_x_continuous(breaks= c(50000,150000),labels = c("", ""))+
  scale_color_manual(values=c(all_read_color,ful_read_color,nfl_read_color))


pdf(file=here("figures", "figure1","figure1B.pdf"),width=4.167,height=2.577,colormodel = "rgb")
p4
dev.off()