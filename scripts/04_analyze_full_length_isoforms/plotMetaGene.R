#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
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
# 
# #---L1-1--
# 
# all_data <-   read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.all_meta.txt"),header=FALSE)
# full_data <-  read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.full_meta.txt"),header=FALSE)
# short_data <- read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.short_meta.txt"),header=FALSE)
# 
# df <- data.frame(x = 1:200000,
#                  all_reads = all_data,
#                  full_reads = full_data,
#                  short_reads = short_data)
# colnames(df) <- c("x","all_reads","full_reads","short_reads")
# 
# p4 <- ggplot(df) +
#   geom_vline(xintercept=50000,size=0.5)+
#   geom_vline(xintercept=150000,size=0.5)+
#   geom_line(aes(x=x, y=all_reads,color = "all reads"),size = 1)+
#   geom_line(aes(x=x, y=full_reads,color = "full-length\nreads"),size = 1)+
#   geom_line(aes(x=x, y=short_reads,color = "non full-length\nreads"),size = 1)+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x= element_blank(),
#         axis.text.y = element_text(size=10,family="Helvetica"),
#         axis.title.y= element_text(size=10,family="Helvetica"),
#         legend.title = element_blank(),
#         legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
#         legend.key.size = unit(1.25, 'lines'),
#         legend.position = "top",
#         legend.justification = c(0.5, 0))+
#   geom_rect(data=data.frame(),aes(xmin=1, xmax=200000),ymin=-2.5e-7 - 3.5e-7- 1e-7,ymax=2.5e-7 - 3.5e-7- 1e-7,fill=gene_model_color,size =0.5,color="black")+
#   geom_rect(data=data.frame(),aes(xmin=50000,xmax=150000),ymin=-3.5e-7 - 3.5e-7 - 2e-7,  ymax=3.5e-7 - 3.5e-7,  fill=gene_model_color,size =0.5,color="black")+
#   geom_text(data=data.frame(x=1),x=25000, y=-7e-7 - 3.5e-7- 3e-7,label = "5'UTR",family="Helvetica",size=3.52734)+
#   geom_text(data=data.frame(x=1),x=100000, y=-7e-7 - 3.5e-7- 3e-7,label = "CDS",family="Helvetica",size=3.52734)+
#   geom_text(data=data.frame(x=1),x=175000, y=-7e-7 - 3.5e-7- 3e-7,label = "3'UTR",family="Helvetica",size=3.52734)+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0))+
#   coord_cartesian(ylim=c(-12e-7- 2e-7,8.25e-6),clip="off")+
#   ylab("Relative Coverage")+
#   scale_x_continuous(breaks= c(50000,150000),labels = c("", ""))+
#   scale_color_manual(values=c(all_read_color,ful_read_color,nfl_read_color))
# 
# 
# pdf(file=here("figures", "figure1","figure1B.pdf"),width=4.167,height=2.577,colormodel = "rgb")
# p4
# dev.off()


# #---L1-1--
# 
# all_data <-   read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.all_meta.txt"),header=FALSE)
# full_data <-  read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.only_stringent_meta.txt"),header=FALSE)
# short_data <- read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.only_stringent_short_meta.txt"),header=FALSE)
# 
# df <- data.frame(x = 1:200000,
#                  all_reads = all_data,
#                  full_reads = full_data,
#                  short_reads = short_data)
# colnames(df) <- c("x","all_reads","full_reads","short_reads")
# 
# p4 <- ggplot(df) +
#   geom_vline(xintercept=50000,size=0.5)+
#   geom_vline(xintercept=150000,size=0.5)+
#   geom_line(aes(x=x, y=all_reads,color = "all reads"),size = 1)+
#   geom_line(aes(x=x, y=short_reads,color = "non-stringent\nreads"),size = 1)+
#   geom_line(aes(x=x, y=full_reads,color = "stringent\nreads"),size = 1)+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x= element_blank(),
#         axis.text.y = element_text(size=10,family="Helvetica"),
#         axis.title.y= element_text(size=10,family="Helvetica"),
#         legend.title = element_blank(),
#         legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
#         legend.key.size = unit(1.25, 'lines'),
#         legend.position = "top",
#         # legend.position = c(0.5,0),
#         legend.justification = c(0.5, 0))+
#   geom_rect(data=data.frame(),aes(xmin=1, xmax=200000),ymin=-2.5e-7 - 3.5e-7- 1e-7,ymax=2.5e-7 - 3.5e-7- 1e-7,fill=gene_model_color,size =0.5,color="black")+
#   geom_rect(data=data.frame(),aes(xmin=50000,xmax=150000),ymin=-3.5e-7 - 3.5e-7 - 2e-7,  ymax=3.5e-7 - 3.5e-7,  fill=gene_model_color,size =0.5,color="black")+
#   geom_text(data=data.frame(x=1),x=25000, y=-7e-7 - 3.5e-7- 3e-7,label = "5'UTR",family="Helvetica",size=3.52734)+
#   geom_text(data=data.frame(x=1),x=100000, y=-7e-7 - 3.5e-7- 3e-7,label = "CDS",family="Helvetica",size=3.52734)+
#   geom_text(data=data.frame(x=1),x=175000, y=-7e-7 - 3.5e-7- 3e-7,label = "3'UTR",family="Helvetica",size=3.52734)+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0))+
#   coord_cartesian(ylim=c(-12e-7- 2e-7,8.25e-6),clip="off")+
#   ylab("Relative Coverage")+
#   scale_x_continuous(breaks= c(50000,150000),labels = c("", ""))+
#   scale_color_manual(values=c(all_read_color,nfl_read_color,ful_read_color))
# 
# 
# pdf(file=here("figures", "figure1","figure1B_v3.pdf"),width=4.167,height=2.577,colormodel = "rgb")
# p4
# dev.off()
# 
# 
# ### 
# 
# all_data <-   read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.all_meta.txt"),header=FALSE)
# full_data <-  read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.short_meta.txt"),header=FALSE)
# short_data <- read.table(here("data","scratch","metagene","ce11_gen_L1_bio1_tech1_Minimap2_k14.only_stringent_short_meta.txt"),header=FALSE)
# 
# df <- data.frame(x = 1:200000,
#                  all_reads = all_data,
#                  full_reads = full_data,
#                  short_reads = short_data)
# colnames(df) <- c("x","all_reads","full_reads","short_reads")
# 
# p4 <- ggplot(df) +
#   geom_vline(xintercept=50000,size=0.5)+
#   geom_vline(xintercept=150000,size=0.5)+
#   geom_line(aes(x=x, y=all_reads,color = "all reads"),size = 1)+
#   geom_line(aes(x=x, y=short_reads,color = "non-stringent\nreads"),size = 1)+
#   geom_line(aes(x=x, y=full_reads,color = "non-full length\nreads"),size = 1)+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.line.x = element_blank(),
#         axis.ticks.x= element_blank(),
#         axis.text.y = element_text(size=10,family="Helvetica"),
#         axis.title.y= element_text(size=10,family="Helvetica"),
#         legend.title = element_blank(),
#         legend.text = element_text(size=10,family="Helvetica",margin=margin(r=8,l=4)),
#         legend.key.size = unit(1.25, 'lines'),
#         legend.position = "top",
#         # legend.position = c(0.5,0),
#         legend.justification = c(0.5, 0))+
#   geom_rect(data=data.frame(),aes(xmin=1, xmax=200000),ymin=-2.5e-7 - 3.5e-7- 1e-7,ymax=2.5e-7 - 3.5e-7- 1e-7,fill=gene_model_color,size =0.5,color="black")+
#   geom_rect(data=data.frame(),aes(xmin=50000,xmax=150000),ymin=-3.5e-7 - 3.5e-7 - 2e-7,  ymax=3.5e-7 - 3.5e-7,  fill=gene_model_color,size =0.5,color="black")+
#   geom_text(data=data.frame(x=1),x=25000, y=-7e-7 - 3.5e-7- 3e-7,label = "5'UTR",family="Helvetica",size=3.52734)+
#   geom_text(data=data.frame(x=1),x=100000, y=-7e-7 - 3.5e-7- 3e-7,label = "CDS",family="Helvetica",size=3.52734)+
#   geom_text(data=data.frame(x=1),x=175000, y=-7e-7 - 3.5e-7- 3e-7,label = "3'UTR",family="Helvetica",size=3.52734)+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0))+
#   coord_cartesian(ylim=c(-12e-7- 2e-7,8.25e-6),clip="off")+
#   ylab("Relative Coverage")+
#   scale_x_continuous(breaks= c(50000,150000),labels = c("", ""))+
#   scale_color_manual(values=c(all_read_color,nfl_read_color,ful_read_color))
# 
# 
# pdf(file=here("figures", "figure1","figure1B_v4.pdf"),width=4.167,height=2.577,colormodel = "rgb")
# p4
# dev.off()

###

all_data <-   read.table(here("data","scratch","metagene","L1.original_alignment_meta.txt"),header=FALSE)
full_data <-  read.table(here("data","scratch","metagene","L1.beta_pass_meta.txt"),header=FALSE)
short_data <- read.table(here("data","scratch","metagene","L1.beta_fail_meta.txt"),header=FALSE)

df <- data.frame(x = 1:200000,
                 all_reads = all_data,
                 full_reads = full_data,
                 short_reads = short_data)
colnames(df) <- c("x","all_reads","full_reads","short_reads")

p4 <- ggplot(df) +
  geom_vline(xintercept=50000,size=0.5)+
  geom_vline(xintercept=150000,size=0.5)+
  geom_line(aes(x=x, y=all_reads,color = "all reads"),size = 1)+
  geom_line(aes(x=x, y=short_reads,color = "non full-length\nreads"),size = 1)+
  geom_line(aes(x=x, y=full_reads,color = "full-length\nreads"),size = 1)+
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
        # legend.position = c(0.5,0),
        legend.justification = c(0.5, 0))+
  geom_rect(data=data.frame(),aes(xmin=1, xmax=200000),ymin=-2.5e-7 - 3.5e-7- 1e-7,ymax=2.5e-7 - 3.5e-7- 1e-7,fill=gene_model_color,size =0.5,color="black")+
  geom_rect(data=data.frame(),aes(xmin=50000,xmax=150000),ymin=-3.5e-7 - 3.5e-7 - 2e-7,  ymax=3.5e-7 - 3.5e-7,  fill=gene_model_color,size =0.5,color="black")+
  geom_text(data=data.frame(x=1),x=25000, y=-7e-7 - 3.5e-7- 3e-7,label = "5'UTR",family="Helvetica",size=3.52734)+
  geom_text(data=data.frame(x=1),x=100000, y=-7e-7 - 3.5e-7- 3e-7,label = "CDS",family="Helvetica",size=3.52734)+
  geom_text(data=data.frame(x=1),x=175000, y=-7e-7 - 3.5e-7- 3e-7,label = "3'UTR",family="Helvetica",size=3.52734)+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=c(-12e-7- 2e-7,8.25e-6),clip="off")+
  ylab("Relative Coverage")+
  scale_x_continuous(breaks= c(50000,150000),labels = c("", ""))+
  scale_color_manual(values=c(all_read_color,ful_read_color,nfl_read_color))


pdf(file=here("figures", "figure1","figure1B.pdf"),width=4.167,height=2.577,colormodel = "rgb")
p4
dev.off()
