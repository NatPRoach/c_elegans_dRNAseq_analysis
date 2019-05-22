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
df <- read.table("results/scratch/countGenesAndIsoforms/full_length_overlap.matrix",header=TRUE,sep='\t')
g1 <- ggplot(df,aes(class))+
  geom_bar(aes(fill = support), position = position_stack(reverse = TRUE),width=0.25,color="black")+
  scale_fill_manual(values=c(full_length_color,inferred_color))+
  scale_x_discrete(labels=c("Our genes" = "genes", "Our isoforms" = "isoforms","WB genes"="WB\ngenes","WB isoforms"="WB\nisoforms"))+
  ylab("Number identified")+
  theme(text = element_text(size = 8,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,family = "Helvetica"),
        axis.text.y = element_text(size=8,family = "Helvetica"),
        legend.text = element_text(size=8,family = "Helvetica"),
        legend.title = element_blank())


pdf(file="figures/figure2/figure2A.pdf",height=2,width=4.5,colormodel="rgb")
print(g1)
dev.off()

# Uncomment to plot venn diagrams for figure 2B
gene_overlap <- read.table("results/scratch/countGenesAndIsoforms/full_length_gene_overlap.matrix",sep="\t",header = FALSE)
colnames(gene_overlap) <- c("supported\ngenes","supported\nWB genes")
e1 <- euler(gene_overlap)
pdf(file="figures/supplementals/sfigure3/sfigure3A.pdf",height=2,width=3,colormodel="rgb")
plot(e1,fills = c(supported_color,supported_WB_color),labels=list(fontsize=8),quantities=list(fontsize=8))
dev.off()

isoform_overlap <- read.table("results/scratch/countGenesAndIsoforms/full_length_isoform_overlap.matrix",sep="\t",header = FALSE)
colnames(isoform_overlap) <- c("supported\nisoforms","supported\nWB isoforms")
e2 <- euler(isoform_overlap)
pdf(file="figures/figure2/figure2B.pdf",height=2,width=3,colormodel="rgb")
plot(e2,fills = c(supported_color,supported_WB_color),labels=list(fontsize=8),quantities=list(fontsize=8))
dev.off()


df <- read.table("results/scratch/countGenesAndIsoforms/staged_gene_and_isoform_count.txt",header=TRUE,sep='\t')
df$dataset <- factor(df$dataset,levels = c("Our genes","Our isoforms"))
df <- df[df$dataset == "Our genes" | df$dataset == "Our isoforms",]
g2 <- ggplot(df,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(genes_color,splice_isoform_color,"#E01E26","#FD7F23"),labels=c("Our genes" = "genes", "Our isoforms"= "splice\nisoforms"))+
  ylab("Number identified")+
  xlab("Stage")+
  theme(text = element_text(size = 8,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,family = "Helvetica"),
        axis.text.y = element_text(size=8,family = "Helvetica"),
        legend.text = element_text(size=8,family = "Helvetica"),
        legend.title = element_blank())

pdf(file="figures/figure2/figure2C.pdf",height=2,width=4.5,colormodel="rgb")
print(g2)
dev.off()

df2 <- read.table("results/scratch/countGenesAndIsoforms/staged_novel_gene_and_isoform_count.txt",header=TRUE,sep='\t')
g3 <- ggplot(df2,aes(x=reorder(stage,x_order),y=counts))+
  geom_bar(position="dodge",stat = "identity",aes(fill = dataset),color="black")+
  scale_x_discrete(labels = c("L1"= "L1","L2"= "L2","L3"= "L3","L4"= "L4","young adult"= "young\nadult","mature adult"= "mature\nadult","male"= "male","all"= "all"))+
  scale_fill_manual(values=c(genes_with_novel_isoforms_color,novel_isoforms_color),labels=c("genes with\nnovel isoforms","novel isoforms"))+
  ylab("Number identified")+
  xlab("Stage")+
  theme(text = element_text(size = 8,family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,family = "Helvetica"),
        axis.text.y = element_text(size=8,family = "Helvetica"),
        legend.text = element_text(size=8,family = "Helvetica"),
        legend.title = element_blank())


pdf(file="figures/figure2/figure2E.pdf",height=2,width=4.5,colormodel="rgb")
print(g3)
dev.off()

