#!/Library/Frameworks/R.framework/Resources/bin/Rscript
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(eulerr)
library(here)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
full_length_color = "#389E34"
inferred_color = "#2579B2"
supported_color = "#A7CEE2"
supported_WB_color = "#B3DE8E"
genes_color = "#A7CEE2"
splice_isoform_color = "#B3DE8E"
genes_with_novel_isoforms_color = "#A7CEE2"
novel_isoforms_color = "#B3DE8E"
}else{
  full_length_color = args[1]
  inferred_color = args[2]
  supported_color = args[3]
  supported_WB_color = args[4]
  genes_color = args[5]
  splice_isoform_color = args[6]
  genes_with_novel_isoforms_color = args[7]
  novel_isoforms_color = args[8]
}


setwd(here())


# df <- read.table("results/scratch/countGenesAndIsoforms/full_length_overlap.matrix",header=TRUE,sep='\t')
# g1 <- ggplot(df,aes(class))+
#   geom_bar(aes(fill = support), position = position_stack(reverse = TRUE),width=0.25,color="black")+
#   scale_fill_manual(values=c(full_length_color,inferred_color))+
#   scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms","WB genes"="WB\ngenes","WB isoforms"="WB\nisoforms"))+
#   ylab("Number identified")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())
# 
# 
# pdf(file="figures/figure2/figure2A.pdf",height=2,width=4.5,colormodel="rgb")
# print(g1)
# dev.off()
# 
# # Uncomment to plot venn diagrams for figure 2B
# gene_overlap <- read.table("results/scratch/countGenesAndIsoforms/full_length_gene_overlap.matrix",sep="\t",header = FALSE)
# colnames(gene_overlap) <- c("supported\ngenes","supported\nWB genes")
# e1 <- euler(gene_overlap)
# pdf(file="figures/supplementals/sfigure3/sfigure3A.pdf",height=2,width=3,colormodel="rgb")
# plot(e1,fills = c(supported_color,supported_WB_color),labels=list(fontsize=8),quantities=list(fontsize=8))
# dev.off()
# 
# isoform_overlap <- read.table("results/scratch/countGenesAndIsoforms/full_length_isoform_overlap.matrix",sep="\t",header = FALSE)
# colnames(isoform_overlap) <- c("supported\nisoforms","supported\nWB isoforms")
# e2 <- euler(isoform_overlap)
# pdf(file="figures/figure2/figure2B.pdf",height=2,width=3,colormodel="rgb")
# plot(e2,fills = c(supported_color,supported_WB_color),labels=list(fontsize=8),quantities=list(fontsize=8))
# dev.off()
# 
# 
# df <- read.table("results/scratch/countGenesAndIsoforms/staged_gene_and_isoform_count.txt",header=TRUE,sep='\t')
# df$dataset <- factor(df$dataset,levels = c("Our genes","Our isoforms"))
# df <- df[df$dataset == "Our genes" | df$dataset == "Our isoforms",]
# g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
#   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice\nisoforms"))+
#   ylab("Number identified")+
#   xlab("Stage")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())
# 
# pdf(file="figures/figure2/figure2C.pdf",height=2,width=4.5,colormodel="rgb")
# print(g2)
# dev.off()
# 
# df2 <- read.table("results/scratch/countGenesAndIsoforms/staged_novel_gene_and_isoform_count.txt",header=TRUE,sep='\t')
# g3 <- ggplot(df2,aes(x=reorder(stage,x_order),y=counts))+
#   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_fill_manual(values=c(genes_with_novel_isoforms_color,novel_isoforms_color),labels=c("genes with\nnovel isoforms","novel isoforms"))+
#   ylab("Number identified")+
#   xlab("Stage")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())
# 
# 
# pdf(file="figures/figure2/figure2E.pdf",height=2,width=4.5,colormodel="rgb")
# print(g3)
# dev.off()


# ### Sensitive 
# 
# df <- read.table("results/scratch/countGenesAndIsoforms/sensitive_full_length_overlap.matrix",header=TRUE,sep='\t')
# # g1 <- ggplot(df,aes(class))+
# #   geom_bar(aes(fill = support), position = position_stack(reverse = TRUE),width=0.25,color="black")+
# #   scale_fill_manual(values=c(full_length_color,inferred_color))+
# #   scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms","WB genes"="WB\ngenes","WB isoforms"="WB\nisoforms"))+
# #   ylab("Number identified")+
# #   theme(text = element_text(size = 8,family = "Helvetica"),
# #         axis.title.x = element_blank(),
# #         axis.text.x = element_text(size=8,family = "Helvetica"),
# #         axis.text.y = element_text(size=8,family = "Helvetica"),
# #         legend.text = element_text(size=8,family = "Helvetica"),
# #         legend.title = element_blank())
# g1 <- ggplot(df,aes(class))+
#   geom_bar(aes(fill = support), position = position_stack(reverse = TRUE),width=0.25,color="black")+
#   scale_fill_manual(values=c(full_length_color,inferred_color))+
#   #scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms"))+
#   scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms","WB genes"="WB\ngenes","WB isoforms"="WB\nisoforms"))+
#   ylab("Number identified")+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0))+
#   theme(text = element_text(size = 10,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=10,family = "Helvetica"),
#         axis.text.y = element_text(size=10,family = "Helvetica"),
#         legend.text = element_text(size=10,family = "Helvetica",margin=margin(r=8,l=4)),
#         legend.title = element_blank(),
#         #legend.position = "top",
#         legend.position = c(0.5,0.90),
#         legend.direction = "horizontal",
#         #legend.justification = c(0.5, 0),
#         legend.justification = c(0.5, 0))
# 
# # pdf(file="figures/figure2/sensitive_figure2A.pdf",height=2,width=4.5,colormodel="rgb")
# pdf(file="figures/figure2/sensitive_figure2A.pdf",height=2.5,width=4.1667,colormodel="rgb")
# 
# print(g1)
# dev.off()
# 
# # Uncomment to plot venn diagrams for figure 2B
# gene_overlap <- read.table("results/scratch/countGenesAndIsoforms/sensitive_full_length_gene_overlap.matrix",sep="\t",header = FALSE)
# colnames(gene_overlap) <- c("supported\ngenes","supported\nWB genes")
# e1 <- euler(gene_overlap)
# pdf(file="figures/supplementals/sfigure3/sensitive_sfigure3A.pdf",height=2.5,width=4.1667,colormodel="rgb")
# plot(e1,fills = c(supported_color,supported_WB_color),labels=list(fontsize=10),quantities=list(fontsize=10))
# dev.off()
# 
# isoform_overlap <- read.table("results/scratch/countGenesAndIsoforms/sensitive_full_length_isoform_overlap.matrix",sep="\t",header = FALSE)
# colnames(isoform_overlap) <- c("supported\nisoforms","supported\nWB isoforms")
# e2 <- euler(isoform_overlap)
# pdf(file="figures/figure2/sensitive_figure2B.pdf",height=2.5,width=4.1667,colormodel="rgb")
# plot(e2,fills = c(supported_color,supported_WB_color),labels=list(fontsize=10),quantities=list(fontsize=10))
# dev.off()
# 
# 
# df <- read.table("results/scratch/countGenesAndIsoforms/sensitive_staged_gene_and_isoform_count.txt",header=TRUE,sep='\t')
# df$dataset <- factor(df$dataset,levels = c("Our genes","Our isoforms"))
# df <- df[df$dataset == "Our genes" | df$dataset == "Our isoforms",]
# # g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
# #   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
# #   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
# #   scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice\nisoforms"))+
# #   ylab("Number identified")+
# #   xlab("Stage")+
# #   theme(text = element_text(size = 8,family = "Helvetica"),
# #         axis.title.x = element_blank(),
# #         axis.text.x = element_text(size=8,family = "Helvetica"),
# #         axis.text.y = element_text(size=8,family = "Helvetica"),
# #         legend.text = element_text(size=8,family = "Helvetica"),
# #         legend.title = element_blank())
# 
# g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
#   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
#   #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all"))+
#   #scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice\nisoforms"))+
#   scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice isoforms"))+
#   ylab("Number identified")+
#   xlab("Stage")+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0))+
#   theme(text = element_text(size = 10,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=10,family = "Helvetica"),
#         axis.text.y = element_text(size=10,family = "Helvetica"),
#         legend.text = element_text(size=10,family = "Helvetica",margin=margin(r=8,l=4)),
#         legend.title = element_blank(),
#         #legend.position = "top",
#         legend.position = c(0.5,0.85),
#         legend.direction = "horizontal",
#         legend.justification = c(0.5, 0))
# 
# # pdf(file="figures/figure2/sensitive_figure2C.pdf",height=2,width=4.5,colormodel="rgb")
# pdf(file="figures/figure2/sensitive_figure2C.pdf",height=2.5,width=4.1667,colormodel="rgb")
# print(g2)
# dev.off()
# 
# df2 <- read.table("results/scratch/countGenesAndIsoforms/sensitive_staged_novel_gene_and_isoform_count.txt",header=TRUE,sep='\t')
# # g3 <- ggplot(df2,aes(x=reorder(stage,x_order),y=counts))+
# #   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
# #   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
# #   scale_fill_manual(values=c(genes_with_novel_isoforms_color,novel_isoforms_color),labels=c("genes with\nnovel isoforms","novel isoforms"))+
# #   ylab("Number identified")+
# #   xlab("Stage")+
# #   theme(text = element_text(size = 8,family = "Helvetica"),
# #         axis.title.x = element_blank(),
# #         axis.text.x = element_text(size=8,family = "Helvetica"),
# #         axis.text.y = element_text(size=8,family = "Helvetica"),
# #         legend.text = element_text(size=8,family = "Helvetica"),
# #         legend.title = element_blank())
# g3 <- ggplot(df2,aes(x=reorder(stage,x_order),y=counts))+
#   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
#   #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all"))+
#   scale_fill_manual(values=c(genes_with_novel_isoforms_color,novel_isoforms_color),labels=c("genes with\nnovel isoforms","novel isoforms"))+
#   #scale_fill_manual(values=c("#2579B2","#389E34"))+
#   ylab("Number identified")+
#   xlab("Stage")+
#   expand_limits(y=0)+
#   scale_y_continuous(expand=c(0,0))+
#   theme(text = element_text(size = 10,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=10,family = "Helvetica"),
#         axis.text.y = element_text(size=10,family = "Helvetica"),
#         legend.text = element_text(size=10,family = "Helvetica",margin=margin(r=8,l=4)),
#         legend.title = element_blank(),
#         #legend.position = "top",
#         legend.position = c(0.5,0.85),
#         legend.direction = "horizontal",
#         legend.justification = c(0.5, 0))
# 
# # pdf(file="figures/figure2/sensitive_figure2E.pdf",height=2,width=4.5,colormodel="rgb")
# pdf(file="figures/figure2/sensitive_figure2E.pdf",height=2.5,width=4.1667,colormodel="rgb")
# 
# print(g3)
# dev.off()


### Stringent 

df <- read.table("results/scratch/countGenesAndIsoforms/stringent_full_length_overlap.matrix",header=TRUE,sep='\t')
# g1 <- ggplot(df,aes(class))+
#   geom_bar(aes(fill = support), position = position_stack(reverse = TRUE),width=0.25,color="black")+
#   scale_fill_manual(values=c(full_length_color,inferred_color))+
#   scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms","WB genes"="WB\ngenes","WB isoforms"="WB\nisoforms"))+
#   ylab("Number identified")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())
g1 <- ggplot(df,aes(class))+
  geom_bar(aes(fill = support), position = position_stack(reverse = TRUE),width=0.25,color="black")+
  scale_fill_manual(values=c(full_length_color,inferred_color))+
  #scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms"))+
  scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms","WB genes"="WB\ngenes","WB isoforms"="WB\nisoforms"))+
  ylab("Number identified")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family = "Helvetica",margin=margin(r=8,l=4)),
        legend.title = element_blank(),
        #legend.position = "top",
        legend.position = c(0.5,0.90),
        legend.direction = "horizontal",
        #legend.justification = c(0.5, 0),
        legend.justification = c(0.5, 0))

# pdf(file="figures/figure2/stringent_figure2A.pdf",height=2,width=4.5,colormodel="rgb")
# pdf(file="figures/figure2/stringent_figure2A.pdf",height=2.5,width=4.1667,colormodel="rgb")
pdf(file="figures/figure2/figure2A.pdf",height=2.5,width=4.1667,colormodel="rgb")

print(g1)
dev.off()

# Uncomment to plot venn diagrams for figure 2B
gene_overlap <- read.table("results/scratch/countGenesAndIsoforms/stringent_full_length_gene_overlap.matrix",sep="\t",header = FALSE)
us_genes   = length(gene_overlap[gene_overlap$V1 == '1' & gene_overlap$V2 == '0',]$V1)
them_genes = length(gene_overlap[gene_overlap$V1 == '0' & gene_overlap$V2 == '1',]$V1)
both_genes = length(gene_overlap[gene_overlap$V1 == '1' & gene_overlap$V2 == '1',]$V1)
isoform_overlap <- read.table("results/scratch/countGenesAndIsoforms/stringent_full_length_isoform_overlap.matrix",sep="\t",header = FALSE)
us_isoforms   = length(isoform_overlap[isoform_overlap$V1 == '1' & isoform_overlap$V2 == '0',]$V1)
them_isoforms = length(isoform_overlap[isoform_overlap$V1 == '0' & isoform_overlap$V2 == '1',]$V1)
both_isoforms = length(isoform_overlap[isoform_overlap$V1 == '1' & isoform_overlap$V2 == '1',]$V1)

df <- data.frame(type = c("genes","genes","genes","isoforms","isoforms","isoforms"),
                 value = c(us_genes,both_genes,them_genes,us_isoforms,both_isoforms,them_isoforms),
                 class = c("This study","Both", "WormBase","This study","Both", "WormBase"),
                 order = c(0,2,1,0,2,1))

g <- ggplot(df, aes(x = type,y=value,fill=factor(class,levels=c("This study","Both","WormBase")),order=order))+
  geom_bar(stat="identity",width=0.25,color="black")+
  # geom_text(aes(label=value), size = 10, position = position_stack(vjust = 0.5) )+
  geom_text(aes(label=value), size = 8/14*5, position = position_stack(vjust = 0.5) )+
  scale_y_continuous(expand=c(0,0))+
  scale_x_discrete(labels=c("genes" = "genes\n ", "isoforms" = "isoforms\n "))+
  # scale_fill_manual(values=c("#80B1D3","#B3DE69","#FFED6F"))+
  scale_fill_manual(breaks = c("This study","WormBase","Both"),values=c("#67A9CF","#B3DE8E","#DBCE85"))+
  ylab("Number with \nfull-length support")+
  
  theme(text = element_text(size = 10,family = "Helvetica"),
         axis.title.x = element_blank(),
         axis.text.x = element_text(size=10,family = "Helvetica"),
         axis.text.y = element_text(size=10,family = "Helvetica"),
         legend.text = element_text(size=10,family = "Helvetica",margin=margin(r=8,l=4)),
         legend.title = element_blank(),
         legend.position = "top",
         # legend.position = c(0.5,0.90),
         legend.direction = "horizontal",
         #legend.justification = c(0.5, 0),
         legend.justification = c(0.5, 0))

# pdf(file="figures/figure2/stringent_figure2B_v2.pdf",height=2.5,width=4.1667,colormodel="rgb")
pdf(file="figures/figure2/figure2B.pdf",height=2.5,width=4.1667,colormodel="rgb")
g
dev.off()

# colnames(gene_overlap) <- c("supported\ngenes","supported\nWB genes")
# us = length(gene_overlap[gene_overlap$`supported
# genes` == '1' & gene_overlap$`supported
# WB genes` == '0',]$`supported
# genes`)


# e1 <- euler(gene_overlap)
# pdf(file="figures/supplementals/sfigure3/stringent_sfigure3A.pdf",height=2.5,width=4.1667,colormodel="rgb")
# plot(e1,fills = c(supported_color,supported_WB_color),labels=list(fontsize=10),quantities=list(fontsize=10))
# dev.off()
# 
# isoform_overlap <- read.table("results/scratch/countGenesAndIsoforms/stringent_full_length_isoform_overlap.matrix",sep="\t",header = FALSE)
# colnames(isoform_overlap) <- c("supported\nisoforms","supported\nWB isoforms")
# e2 <- euler(isoform_overlap)
# pdf(file="figures/figure2/stringent_figure2B.pdf",height=2.5,width=4.1667,colormodel="rgb")
# plot(e2,fills = c(supported_color,supported_WB_color),labels=list(fontsize=10),quantities=list(fontsize=10))
# dev.off()





df <- read.table("results/scratch/countGenesAndIsoforms/stringent_staged_gene_and_isoform_count.txt",header=TRUE,sep='\t')
df$dataset <- factor(df$dataset,levels = c("Our genes","Our isoforms"))
df <- df[df$dataset == "Our genes" | df$dataset == "Our isoforms",]
# g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
#   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice\nisoforms"))+
#   ylab("Number identified")+
#   xlab("Stage")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())

g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all"))+
  #scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice\nisoforms"))+
  scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice isoforms"))+
  ylab("Number identified")+
  xlab("Stage")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family = "Helvetica",margin=margin(r=8,l=4)),
        legend.title = element_blank(),
        #legend.position = "top",
        legend.position = c(0.5,0.85),
        legend.direction = "horizontal",
        legend.justification = c(0.5, 0))

# pdf(file="figures/figure2/stringent_figure2C.pdf",height=2,width=4.5,colormodel="rgb")
# pdf(file="figures/figure2/stringent_figure2C.pdf",height=2.5,width=4.1667,colormodel="rgb")
pdf(file="figures/figure2/figure2C.pdf",height=2.5,width=4.1667,colormodel="rgb")
print(g2)
dev.off()

df2 <- read.table("results/scratch/countGenesAndIsoforms/stringent_staged_novel_gene_and_isoform_count.txt",header=TRUE,sep='\t')
# g3 <- ggplot(df2,aes(x=reorder(stage,x_order),y=counts))+
#   geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
#   scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
#   scale_fill_manual(values=c(genes_with_novel_isoforms_color,novel_isoforms_color),labels=c("genes with\nnovel isoforms","novel isoforms"))+
#   ylab("Number identified")+
#   xlab("Stage")+
#   theme(text = element_text(size = 8,family = "Helvetica"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=8,family = "Helvetica"),
#         axis.text.y = element_text(size=8,family = "Helvetica"),
#         legend.text = element_text(size=8,family = "Helvetica"),
#         legend.title = element_blank())
g3 <- ggplot(df2,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  #scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "yAd","mature adult"= "mAd","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(genes_with_novel_isoforms_color,novel_isoforms_color),labels=c("genes with\nnovel isoforms","novel isoforms"))+
  #scale_fill_manual(values=c("#2579B2","#389E34"))+
  ylab("Number identified")+
  xlab("Stage")+
  expand_limits(y=0)+
  scale_y_continuous(expand=c(0,0))+
  theme(text = element_text(size = 10,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=10,family = "Helvetica"),
        axis.text.y = element_text(size=10,family = "Helvetica"),
        legend.text = element_text(size=10,family = "Helvetica",margin=margin(r=8,l=4)),
        legend.title = element_blank(),
        #legend.position = "top",
        legend.position = c(0.5,0.85),
        legend.direction = "horizontal",
        legend.justification = c(0.5, 0))

# pdf(file="figures/figure2/stringent_figure2E.pdf",height=2,width=4.5,colormodel="rgb")
# pdf(file="figures/figure2/stringent_figure2E.pdf", height=2.5,width=4.1667,colormodel="rgb")
pdf(file="figures/figure2/figure2E.pdf", height=2.5,width=4.1667,colormodel="rgb")
print(g3)
dev.off()
