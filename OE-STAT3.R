
#202408903-OE-ST3
rm(list = ls())
setwd("G:\\IB\\202408-RNA-seq\\R\\OE-STAT3")

# 检查输出文件夹是否存在
folder_path <- "./out"
if (!dir.exists(folder_path)) {
  # 如果文件夹不存在，则创建它
  dir.create(folder_path)
  print(paste("Folder created at", folder_path))
} else {
  print(paste("Folder already exists at", folder_path))
}

#0.读取表达矩阵数据----
# setwd("G:\\IB\\202408-RNA-seq\\R\\OE-STAT4")
exp<-readRDS("G:\\IB\\202408-RNA-seq\\R\\OE-STAT4/out/exp_all_fpkm.rdata")
#提出STAT3相关数据
exprSet<-log2(exp[,c(4:6,13:15)] + 1)   #log2  表达矩阵
#分组信息
group <-data.frame(name=colnames(exprSet),#样品名称
                   group=rep(c('NC', 'OE'), each = 3, length.out = 6) ) #分组信息
contra <-c('OE-NC')#对比信息
adjPval <- 0.05         # 矫正p值-p阈值
aflogFC <- 0.5             # logFC  阈值


#1.0 差异分析----

if(T){
  ##差异limma
  library(limma)
  # 构建分组矩阵--design 
  design <- model.matrix(~0+factor(group$group))
  colnames(design) <- levels(factor(group$group))
  rownames(design) <- colnames(exprSet)
  contrast.matrix <- makeContrasts(contra,levels = design)
  
  fit <- lmFit(exprSet,design) #非线性最小二乘法
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
  DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
  DEG <- na.omit(DEG)
  DEG$Symbol <-rownames(DEG)
}

# 输出差异分析表格
write.table(DEG, file = "./out/all_DEG_limma.txt",sep = "\t", row.names = F, col.names = T, quote = F)

#筛选差异表达基因 adjPval差
DEG$Group <- ifelse((DEG$P.Value <adjPval & abs(DEG$logFC)>aflogFC), ifelse(DEG$logFC>aflogFC,"Up","Down"), "Stable")
print(table(DEG$Group))
write.csv(DEG,"./out/all_DEG_limma.csv")#1 差异太少？DOWN 2 UP 3;0.5 down 14, up 25
#volcano----
gene_volcano <- c("STAT3","STAT4","SLC47A1","SLC7A11","GPX4","CD274")#DEG[DEG$Group != "Stable",]$Symbol#rownames(DEG[abs(DEG$logFC) >1.5,])#火山图展示基因 
#兴趣基因
library(ggplot2);library(ggrepel )
p1 <- ggplot(DEG,aes(logFC, -log10(P.Value)))+ #adj.P.Val
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "skyblue")+
  #geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-aflogFC,aflogFC), linetype = "dashed", color = "skyblue")+
  geom_point(aes(size = AveExpr , #-(adj.P.Val)
                 color = AveExpr , #-log10(adj.P.Val)
                 shape=Group
  ))+
  #scale_color_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_colour_viridis_c(direction = -1,option = "D",alpha = 0.7)+ #默认D direction = -1,G H
  scale_size_continuous(range = c(0,1))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.justification = c(0,1))+
  # 设置图例
  guides(col = 
           guide_colorbar(title = "AveExpr" ,#"-Log10 adj.Pval",
                          ticks.colour = NA,
                          #reverse = T,
                          title.vjust = 0.8,
                          barheight = 8,
                          barwidth = 1),
         size = "none") +
  # 添加标签：
  geom_text_repel(data = DEG[DEG$Symbol %in% gene_volcano,],
                  max.overlaps =30,#length(gene_volcano),# getOption("ggrepel.max.overlaps", default = 20),
                  # 这里的filter很关键，筛选你想要标记的基因
                  aes(label =Symbol),
                  size = 2.5, 
                  #box.padding = 0.5, #字到点的距离 
                  #point.padding = 0.8, #字到点的距离，点周围的空白宽度 
                  min.segment.length = 0.1, 
                  segment.size=0.3, 
                  color = 'black') +
  labs(x="Log2(FC)",y="-Log10(pvalue)") #adj.Pval

p1
ggsave("./out/volcano-D-1.pdf",family="serif",width = 5,height = 4)

ggplot(DEG,aes(logFC, -log10(P.Value)))+ #adj.P.Val
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",linewidth=0.1, color = "black")+
  #geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-aflogFC,aflogFC), linetype = "dashed",linewidth=0.1, color = "black")+
  geom_point(aes(size = abs(logFC)*100,#AveExpr , #-(adj.P.Val)
                 color =Group, #-log10(P.Value),# , #-log10(adj.P.Val)
                 shape=Group
  ))+
  #scale_color_gradientn(values = seq(0,1,0.2), colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  #scale_colour_viridis_c(direction = -1,option = "G",alpha = 0.7)+ #默认D direction = -1,G H
  scale_color_manual(values = c("skyblue","grey","pink"))+ #RColorBrewer::brewer.pal(9,"Paired")
  scale_size_continuous(range = c(0,1))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.justification = c(0,1))+
  # 设置图例
  guides(col = 
           guide_colorbar(title = "-log10(pvalue)" ,#"-Log10 adj.Pval",
                          ticks.colour = NA,
                          #reverse = T,
                          title.vjust = 0.8,
                          barheight = 8,
                          barwidth = 1),
         size = "none") +
  # 添加标签：
  geom_text_repel(data = DEG[DEG$Symbol %in% gene_volcano,],
                  max.overlaps =30,#length(gene_volcano),# getOption("ggrepel.max.overlaps", default = 20),
                  # 这里的filter很关键，筛选你想要标记的基因
                  aes(label =Symbol),
                  size = 2, 
                  #box.padding = 0.5, #字到点的距离 
                  #point.padding = 0.8, #字到点的距离，点周围的空白宽度 
                  min.segment.length = 0.1, 
                  segment.size=0.2, 
                  color = 'black') +
  labs(x="Log2(FC)",y="-Log10(pvalue)") #adj.Pval
ggsave("./out/volcano-2-1.pdf",family="serif",width = 5,height = 4)
##rank----
#https://mp.weixin.qq.com/s/Go1_K6MtONJyhkpbmm69Bg

if(T){
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(ggfun)
}
#DEG <-read.csv("./out/all_DEG_limma.csv",row.names = 1)
deg_result <- DEG %>% #read_delim(file = "DEG_result.txt", col_names = T, delim = "\t") %>%
  dplyr::mutate(rank = -1 * rank(logFC, ties.method= "max")) #

deg_result_top_5 <- deg_result %>%
  dplyr::filter(!is.na(Symbol)) %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::slice_head(n = 5)

deg_result_tail_5 <- deg_result %>%
  dplyr::filter(!is.na(Symbol)) %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::slice_tail(n = 5)
#plot
deg_result %>%
  ggplot() + 
  geom_point(aes(x = rank, y = logFC, color = P.Value, size = abs(logFC))) + 
  scale_x_continuous(breaks = c(-15000, -10000, -5000, 0),
                     labels = c(0, 5000, 10000, 15000)) + 
  scale_color_gradient2(low = "#fb8072",high = "#80b1d3",mid = "#ffffff", midpoint = 0.5, name = "pvalue") + 
  geom_hline(yintercept = 0) + 
  # geom_vline(xintercept = -7500) + 
  geom_text_repel(data = deg_result_top_5,
                  aes(x = rank + 10, y = logFC, label = Symbol),
                  box.padding = 0.5,
                  nudge_x = 10,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "left") + 
  geom_text_repel(data = deg_result_tail_5,
                  aes(x = rank + 10, y = logFC, label = Symbol),
                  box.padding = 0.5,
                  nudge_x = 10,
                  nudge_y = -0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "right") + 
  scale_size(range = c(1,8), name = "log2(Fold Change)") + 
  labs(x = "Rank of differentially expressed genes",
       y = "log2(FoldChange)") + 
  theme_bw() + 
  theme(
    panel.border = element_rect(linewidth = 1.25),
    #legend.background = element_roundrect(color = "#808080",linetype = 1),#边框
    axis.text = element_text(color = "#000000", size = 12.5),
    axis.title = element_text(color = "#000000", size = 15)
  )
#显著性展示？

ggsave(filename = "./out/rank.pdf",
       height = 7,
       width = 9)
##兴趣基因热图heatmap----

library(pheatmap);library(RColorBrewer)
#pheatmap(exprSet)#太多基因，卡
#差异基因 DEG[DEG$Group != "Stable",]$Symbol
gene_volcano =DEG[DEG$Group != "Stable",]$Symbol
hallmark_p <-pheatmap(exprSet[rownames(exprSet) %in% gene_volcano, ], #gene_volcano
                      #exprSet[rownames(exprSet) %in% c(DEG[DEG$Group != "Stable",]$Symbol) , ]##差异基因 
                      scale="column",#row
                      cluster_rows = T,cluster_cols  = F, 
                      fontsize_row = 8,fontsize_col = 10, #5 10
                      show_colnames = T,show_rownames = T,
                      border_color = "white",#cellwidth=10,cellheight = 10,
                      gaps_col= 3, #table(exp_marker1$response)
                      color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(8))# ,#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
                      #annotation_col=annotation_col,annotation_row = anno_row
)
hallmark_p

save_pheatmap_pdf <- function(x, filename, width, height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.02)#修改聚类树线条宽度
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数
save_pheatmap_pdf(hallmark_p, "./out/all_pheatmap-DEG.pdf",3,5) #allgene
#save_pheatmap_pdf(hallmark_p, "./out/DEGheatmap.pdf",4.5,17) #allDEGgene

##某个基因T检验boxplot----
gene <- "CD274"

if(T){
  df_gene <-exprSet[rownames(exprSet) %in% gene,]
  df_gene <-as.data.frame(t(df_gene))
  df_gene <-as.data.frame(scale(df_gene) )  #scale样本
  df_gene$Group <- group$group #c(rep("Control",3),rep("Drug",3))  #分组信息！！！
  colnames(df_gene)[1] <-"Gene"
  library(ggplot2);library(ggpubr)
  ggplot(df_gene,#
         aes(Group,Gene,color=Group))+  #age_group
    geom_boxplot(alpha=1,width=0.45,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    geom_point(alpha=0.7,size=2,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    scale_color_manual(values =  c("#EDA065","#66CCFF","#7EC7A7"))+ #blue c("#0066CC","#66CCFF","#9DC3E6")
    theme_bw() +labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
    theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
          panel.grid = element_blank(),legend.key = element_blank(),legend.position="none" ) +
    # stat_compare_means(#aes(group = Group) ,
    # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
    # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
    # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
    # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
    stat_compare_means(method = "wilcox.test",size=5,show.legend= F,label = "p.format",label.x =0.75,label.y =max(df_gene$Gene)-0.06)#+ #,label.y =max(df_gene$Gene)-0.55,label.x =1.5,label.y = 1.3  0.8 ,label.y = 4 #非参数检验kruskal.test not anova
  ylim(NA,max(df_gene$Gene)+0.5) #wilcox.test  t.test
  ggsave(paste0("./out/boxplot_Ttest_",gene,".pdf"),width = 2,height = 3,family="serif")
  
}
#ggsave(paste0("./out/boxplot_",gene,"_Ttest.pdf"),width = 2,height = 3,family="serif")

##批量基因boxplot----
gene_more <-unique(c(gene_volcano,c("STAT3","STAT4","SLC47A1","SLC7A11","GPX4","CD274"))) #输入基因集
for (i in gene_more){
  print(i)
  df_gene <-exprSet[rownames(exprSet) %in% i,]
  df_gene <-as.data.frame(t(df_gene))
  df_gene <-as.data.frame(scale(df_gene) )  #scale样本
  df_gene$Group <- group$group #c(rep("Control",3),rep("Drug",3))  #分组信息！！！
  colnames(df_gene)[1] <-"Gene"
  library(ggplot2);library(ggpubr)
  ggplot(df_gene,#
         aes(Group,Gene,color=Group))+  #age_group
    geom_boxplot(alpha=1,width=0.45,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    geom_point(alpha=0.7,size=2,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    scale_color_manual(values =  c("#EDA065","#66CCFF","#7EC7A7"))+ #blue c("#0066CC","#66CCFF","#9DC3E6")
    theme_bw() +labs(color="",title = i,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
    theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
          panel.grid = element_blank(),legend.key = element_blank(),legend.position="none" ) +
    # stat_compare_means(#aes(group = Group) ,
    # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
    # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
    # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
    # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
    stat_compare_means(method = "t.test",size=5,show.legend= F,label = "p.format",label.x =0.75,label.y =max(df_gene$Gene)-0.06)#+ #,label.y =max(df_gene$Gene)-0.55,label.x =1.5,label.y = 1.3  0.8 ,label.y = 4 #非参数检验kruskal.test not anova
  ylim(NA,max(df_gene$Gene)+0.5) #wilcox.test更加严格  t.test
  ggsave(paste0("./out/boxplot_Ttest_",i,".pdf"),width = 2,height = 3,family="serif")
  
  
}
##gene_volcano top 热图+显著性----
#T检验显著性表格
pvalue_interest <-data.frame() #T.test
exp_int<-exprSet[rownames(exprSet) %in% gene_volcano,] #gene_volcano gene_more兴趣基因集
for (i in 1:nrow(exp_int) ){ 
  print(i)
  pwilcox <- t.test(exp_int[i,4:6],exp_int[i,1:3]) #样本数
  fc <- mean(as.numeric( exp_int[i,4:6] ) )/mean( as.numeric(exp_int[i,1:3]) )#drug/control  手动输入？ 对照组在前！
  pvalue_interest<-rbind(pvalue_interest,data.frame(pvalue=pwilcox$p.value,Fold=fc,Symbol=rownames(exp_int)[i] ) )
}
pvalue_interest$Sign. <- ifelse(pvalue_interest$pvalue>= 0.05,"ns",
                                ifelse(pvalue_interest$pvalue < 0.001,"***",
                                       ifelse(pvalue_interest$pvalue < 0.01,"**","*")
                                ) 
)
pvalue_interest$FC <-ifelse(pvalue_interest$Fold <=0.5,"Decrease",
                            ifelse(pvalue_interest$Fold >=1,"Increase","Unchange") )
rownames(pvalue_interest) <-pvalue_interest$Symbol
#图注
anno_row_p <- data.frame(
  Sign.=pvalue_interest$Sign.,FoldChange=pvalue_interest$FC) #, size = 20, replace = TRUE  adj.P.Val
rownames(anno_row_p) <- rownames(pvalue_interest)
ann_colors_p=list(FoldChange=c('Unchange'='#CCECFF','Decrease'="skyblue",'Increase'="#4a6fe3"),  #'Unchange'='#CAB2D6',
                  Sign.=c("ns"="grey","*"="pink",'**'="#D33F6A","***"="red" ) )
col<-data.frame(group$group);rownames(col)=group$name;colnames(col)="Group"
int_p <-pheatmap(exp_int,
                 scale="column",#row
                 cluster_rows = T,cluster_cols  = F, 
                 fontsize_row = 8,fontsize_col = 10, #5 10
                 show_colnames = T,show_rownames = T,
                 border_color = "white",#cellwidth=10,cellheight = 10,
                 gaps_col= 3, #table(exp_marker1$response)
                 color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)),#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
                 annotation_col=col,annotation_names_row = F,annotation_row = anno_row_p,anno_names_row = F,annotation_colors = ann_colors_p
)
int_p
save_pheatmap_pdf(int_p, "./out/int_pheatmap-1.pdf",5,5) #allDEGgene

#2.0 富集分析----
#使用差异基因
input_gene <-DEG[DEG$Group != "Stable",]$Symbol #39 #DEG[DEG$Group==gene_type,]$Symbol #输入上下调基因
gene_type<-"ALL"
#data文件夹
enrichment <-function(gene_type,input_gene,expr,group,contra){
  
  #1.0 GO----
  library(clusterProfiler)
  trans = bitr(input_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  library(dplyr)
  eego <- enrichGO(gene          = trans$ENTREZID,
                   #universe      = names(geneList),
                   OrgDb         = "org.Hs.eg.db",
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   #pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
  #plot
  df_go <- data.frame(eego) %>%
    group_by(ONTOLOGY) %>%
    slice_head(n = 10) %>% #前10个,前5个？
    arrange(desc(pvalue))
  ratio <- lapply(df_go$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
  df_go$ratio <- ratio
  df_go$Description <- factor(df_go$Description,levels = df_go$Description)
  
  #plot
  my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  library(ggplot2)
  ggplot(df_go) +
    ggforce::geom_link(aes(x = 0,y = Description,
                           xend = -log10(pvalue),yend = Description,
                           alpha = after_stat(index),#stat(index),
                           color = ONTOLOGY,
                           size = 10),#after_stat(index)),#流星拖尾效果
                       n = 500,
                       #color = "#FF0033",
                       show.legend = F) +
    geom_point(aes(x = -log10(pvalue),y = Description),
               color = "black",
               fill = "white",size = 4,shape = 21) +
    geom_line(aes(x = ratio*100,y = Description,group = 1),
              orientation = "y",linewidth = 1,color = "#FFCC00") + #线条-比例 #FFCC00 #E59CC4
    scale_x_continuous(sec.axis = sec_axis(~./100,
                                           labels = scales::label_percent(),
                                           name = "Percent of geneRatio")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold.italic"),
          axis.text = element_text(color = "black")) +
    ylab("") + xlab("-log10 Pvalue") +
    facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) + #3
    scale_color_brewer(palette = "Set2")#OrRd Set1 Paired
  #scale_color_manual(values = c('#CCE0F5', '#CCC9E6', '#625D9E') )
  
  folder_path <- "./out" # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  
  ggsave(paste0("./out/DEG_",gene_type,"_GO.pdf"),family="serif",width = 4,height = 7)
  
  #2.0 KEGG----
  
  kegg <- enrichKEGG(
    gene          = trans$ENTREZID,
    keyType     = "kegg",
    organism   = 'hsa',
    pvalueCutoff      = 0.05,
    pAdjustMethod     = "BH",
    qvalueCutoff  = 0.05
  )
  library(DOSE)
  kegg<-setReadable(kegg, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")
  
  #2.1 enrichplot
  logFC <- DEG[ DEG$Symbol %in% input_gene,]$logFC
  names(logFC) <- DEG[ DEG$Symbol %in% input_gene,]$Symbol  
  
  library(enrichplot)
  pdf(paste0("./out/DEG_",gene_type,"_KEGG_cnetplot.pdf"),family="serif",width = 6,height = 4)
  enrichplot::cnetplot(kegg,showCategory = 10,#foldChange = logFC,
                       circular=F,colorEdge = TRUE)
  enrichplot::cnetplot(kegg,showCategory = 10,foldChange = logFC,
                       circular=F,colorEdge = TRUE)
  # cnetplot(kegg, showCategory = 10, #选择top10的pathway ，这里也可以用包含pathway名称的向量           
  #          color.params = list(foldChange = logFC, #用上面的logFC值标注基因的颜色
  #                                   edge = T)) #显示pathway的标注
  heatplot(kegg, foldChange=logFC,
           #showCategory = 8,
           symbol = "dot", #"rect" 矩形
           pvalue = NULL,
           label_format = 30)
  dev.off()
  
  
  ego <- kegg@result
  ego <- ego [order(ego$Count,decreasing = T),]
  ego <- ego[ego$pvalue <=0.05,] #P > 0.05
  ego <- ego[,c(2,5,8,9)]
  library(ggplot2);library(ggsci);library("scales")
  ego$Description <- factor(ego$Description,levels = c(ego$Description))
  
  #2.2 barplot
  mytheme<- theme(axis.title = element_text(size = 13),
                  axis.text = element_text(size = 11),
                  plot.title = element_text(size = 14,
                                            hjust= 0.5,
                                            face= "bold"),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 11))
  
  p<- ggplot(data = ego,
             aes(x = Count, y = Description, fill = -log10(pvalue)) )+
    scale_fill_viridis_c(alpha = 0.7)+
    #scale_fill_distiller(palette = "RdPu",direction = -1) + #RdPu YlOrRd PiYG
    ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 
    geom_bar(stat = "identity", width = 0.8) +
    labs(x = "Count", y = "", title = "KEGG") + 
    theme_bw()+ mytheme
  p
  ggsave(paste0("./out/DEG_",gene_type,"_KEGG_bar.pdf"),family="serif",width = 7,height = 4)
  #2.2 dotplot
  p<- ggplot(data = ego,
             aes(x = Count, y = Description) )+ #, fill = -log10(pvalue)
    #scale_fill_distiller(palette = "RdPu",direction = -1) + #RdPu YlOrRd PiYG
    ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 
    #geom_bar(stat = "identity", width = 0.8) +
    geom_point(aes(size=Count,color=-log10(pvalue) ))+
    labs(x = "Count", y = "", title = "KEGG") + 
    theme_bw()+ mytheme+
    scale_color_viridis_c(alpha = 0.7)
  p
  ggsave(paste0("./out/DEG_",gene_type,"_KEGG_dot.pdf"),family="serif",width = 7,height = 4)
  
  #2.4 barplot+基因
  #https://mp.weixin.qq.com/s/vZXxoXSOX6Y7bnIuTHaBzg
  # 先自定义主题：
  mytheme1 <- theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
    axis.ticks.length.y = unit(0,"cm"),
    plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
    #legend.title = element_text(size = 13),
    #legend.text = element_text(size = 11),
    legend.position = "none",
    plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
  )
  library(RColorBrewer)
  p <- ggplot(data = ego, aes(x = -log10(pvalue), y = Description, fill =Description )) +  #, fill = Cluster
    #scale_fill_manual(values =c('#6bb9d2', '#d55640')) +
    #scale_fill_manual(values =my36colors) +
    scale_fill_manual(values =rev( colorRampPalette(brewer.pal(11, "Paired"))(length(ego$Description)) )  ) +
    geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
    scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
    #scale_y_continuous(expand = c(0,-1)) +
    labs(x = "-log10(pvalue)", y = "Pathway", title = "KEGG") +
    # x = 0.61 用数值向量控制文本标签起始位置
    geom_text(size=3.8, aes(x = 0.05, label = Description),hjust = 0) + # hjust = 0,左对齐
    #geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust =2.5, color="skyblue" ) + # hjust = 0,左对齐
    geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust =-1, color="skyblue" ) + # hjust = 0,左对齐
    theme_classic() + 
    mytheme1 
  #ylim(min(ego$Description),max(ego$Description))
  #scale_y_discrete(labels = \(x) sub(pattern = "_", replacement = " ", x, fixed = TRUE)) 
  p
  ggsave(paste0("./out/DEG_",gene_type,"_KEGG_bar_genes.pdf"),width = 5,height = 6.5)#7.5
  
  #3.0 GSEA----
  
  #输入所有基因！！！
  trans_all = bitr(DEG$Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  colnames(trans_all)[1]<-"Symbol"
  allDEG <-merge(DEG,trans_all,by="Symbol")
  
  genelist <-  allDEG[allDEG$ENTREZID,]$logFC#所有基因-18866
  names(genelist) = allDEG$ENTREZID
  genelist = sort(genelist, decreasing = TRUE)#排序
  any(genelist == 0)
  genelist<-genelist[genelist>0] #Error in if (abs(max.ES) > abs(min.ES)) { :    missing value where TRUE/FALSE needed
  kkgsea <- gseKEGG(geneList     = genelist,
                    organism     = 'hsa', 
                    minGSSize    = 1,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,#use_internal_data=T,
                    pAdjustMethod = "none" ) #进行gseKEGG富集分析 
  kkgsea=setReadable(kkgsea,keyType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
  kkgsea_res <- kkgsea@result 
  kkgsea_res_p <- kkgsea_res[kkgsea_res$p.adjust <0.05,]#
  write.csv(kkgsea_res, "./out/GSEA_all_gene.csv") #paste0("./out/gsea_",gene_type,".csv")
  
  library(enrichplot)
  #gseaplot2(kkgsea,geneSetID ="hsa00310",color = "red" ) #
  pdf("./out/enrichplot_top10.pdf")
  gseaplot2(
    kkgsea, #gseaResult object，即GSEA结果
    geneSetID = c(1:10),#"hsa00310",#1,#"hsa05168",#富集的ID编号
    title = "", #标题
    color = "red",#GSEA线条颜色
    base_size = 11,#基础字体大小
    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
    pvalue_table = T, #是否添加 pvalue table
    ES_geom = "line" )#running enrichment score用先还是用点ES_geom = "dot"
  dev.off()
  #丑
  
  #排序
  
  ##3.1 柱状图barplot----
  kkgsea_res_p <- kkgsea_res_p %>% dplyr::arrange(desc(NES))
  kkgsea_res_p$Description <- factor(kkgsea_res_p$Description,levels = kkgsea_res_p$Description)
  
  ggplot(data = kkgsea_res_p, aes(x = NES, y = Description, fill =pvalue )) +  #Description, fill = Cluster
    #scale_fill_viridis_c(alpha = 0.7,option = "G")+ #,direction = -1
    #scale_fill_manual(values =rev( colorRampPalette(brewer.pal(11, "Paired"))(length(kkgsea_res_p$Description)) )  ) +
    scale_fill_distiller(palette = "RdPu",direction = -1) + #RdPu YlOrRd PiYG
    geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
    scale_x_continuous(expand = c(0,0),limits = c(min(kkgsea_res_p$NES)-0.5,max(kkgsea_res_p$NES)+0.5)) + # 调整柱子底部与y轴紧贴
    #scale_y_continuous(expand = c(0,-1)) +
    labs(x = "NES", y = "", title = "GSEA") +
    # x = 0.61 用数值向量控制文本标签起始位置
    #geom_text(size=3.8, aes(x = 0.05, label = Description),hjust = 0) + # hjust = 0,左对齐
    #geom_text(size=3, aes(x = NES,y= Description, label = round(p.adjust,2)), hjust = 0, vjust =2.5, color="navy" ) + # hjust = 0,左对齐
    #geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust =-1, color="skyblue" ) + # hjust = 0,左对齐
    theme_classic() #+ 
  #theme(legend.position = "none")
  #mytheme1 
  ggsave("./out/GSEA_all_barplot.pdf",width = 6.5,height = 4)
  
  ##3.2 GSEAvis----
  #3.2.1 gseaNb
  library(GseaVis)
  #批量输出所有显著性KEGG
  pdf(paste0("./out/GSEA_all_GSEAVIS.pdf"),width = 4,height = 4)
  for (i in kkgsea_res_p$ID){
    gsea_name <-as.character(
      kkgsea_res_p[kkgsea_res_p$ID ==i,]$Description
    )
    
    gsea_plot<-gseaNb(object = kkgsea,
                      geneSetID = i, termWidth = 35,#hsa00310 hsa00480
                      #newGsea = F,addPoint = F,
                      #rmSegment = T, #移除红色标记线
                      #rmHt = T,#移除热图
                      addPval = T,#Add NES and Pvalue
                      pvalX = 0.45,pvalY = 0.8, #调整标签位置和颜色
                      pCol = 'red',pHjust = 0,subPlot = 2)
    print(gsea_plot)
    #ggsave(paste0("./out/GSEA_",gsea_name,".pdf"),width = 4,height = 4)
  }
  dev.off()
  
  # gseaNb(object = kkgsea,
  #        geneSetID = c('hsa00310'), #hsa00310 hsa00480
  #        #newGsea = F,addPoint = F,
  #        #rmSegment = T, #移除红色标记线
  #        #rmHt = T,#移除热图
  #        addPval = T,#Add NES and Pvalue
  #        pvalX = 0.65,pvalY = 0.8, #调整标签位置和颜色
  #        pCol = 'red',pHjust = 0,subPlot = 2)
  
  # pdf("./out/KEGG_Lysine_Glu.pdf",width = 4.5,height = 4.5)
  # gseaNb(object = kkgsea,
  #        geneSetID = c("hsa00310", "hsa00480"),
  #        subPlot = 2,
  #        termWidth = 35,
  #        legend.position = c(0.7,0.85),
  #        addPval = T,#Add NES and Pvalue
  #        pvalX = 0.97,pvalY = 0.7, #调整标签位置和颜色
  #        pCol = 'red',pHjust = 0
  #        )
  # dev.off()
  
  #3.2.2 gseaNb-pheatmap 热图
  
  # gsea_plot<-gseaNb(object = kkgsea,
  #                   geneSetID = i, #hsa00310 hsa00480
  #                   #newGsea = F,addPoint = F,
  #                   #rmSegment = T, #移除红色标记线
  #                   #rmHt = T,#移除热图
  #                   addPval = T,#Add NES and Pvalue
  #                   pvalX = 0.65,pvalY = 0.8, #调整标签位置和颜色
  #                   pCol = 'red',pHjust = 0,subPlot = 2)
  #准备表达矩阵
  #expr <-readRDS("./out/exp_all_fpkm.rdata")[,7:12]
  expr$gene_name <-rownames(expr)
  expr<- expr %>% dplyr::select('gene_name',colnames(expr)[1:dim(expr)[2]-1],everything())
  ##表达量文件第一列必须是基因名
  
  
  pdf(paste0("./out/GSEA_all_GSEAVIS_pheatmap.pdf"),width = 6.5,height = 5.2)
  for (i in kkgsea_res_p$ID){
    gsea_name <-as.character(
      kkgsea_res_p[kkgsea_res_p$ID ==i,]$Description
    )
    
    gsea_plot<-gseaNb(object = kkgsea,newGsea = T,
                      geneSetID = i, #gse@result[["ID"]][1]
                      add.geneExpHt = T,
                      exp = expr,
                      exp.col = c('skyblue','white','pink'),#修改热图颜色
                      ght.geneText.size = 10, #调整热图基因标签大小
                      ght.relHight = 0.3, #修改热图相对高度
                      addPval = T,#Add NES and Pvalue
                      pvalX = 0.6,pvalY = 0.7, #调整标签位置和颜色
                      pCol = 'red',pHjust = 0,
                      subPlot = 2  #只保留上2部分
    )
    print(gsea_plot)
    #ggsave(paste0("./out/GSEA_",gsea_name,".pdf"),width = 4,height = 4)
  }
  dev.off()
  
  #4.0 HALLMARK-GSEA----
  library(msigdbr) 
  hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")#下载基因集
  table(hallmark_genesets$gs_name)
  
  
  fgsea_sets<- hallmark_genesets %>%             #切割分组构建通路基因集list
    split(x = .$gene_symbol,             #分组对象
          f = .$gs_name)                 #切割依据
  
  
  
  genelist_1 <-  allDEG$logFC#所有基因-18866
  names(genelist_1) <- allDEG$Symbol
  genelist_1 = sort(genelist_1, decreasing = TRUE)#排序
  
  library(fgsea)
  fgseaRes<- fgsea(fgsea_sets,            #通路基因集
                   stats = genelist_1 )    #差异变化倍数 输入symbol!
  fgseaRes$pathway <-gsub("HALLMARK_","",fgseaRes$pathway)
  fp=ggplot(fgseaRes %>% 
              filter(pval < 0.05) %>%    #过滤筛选
              head(n= 20), 
            aes(reorder(pathway, NES), NES)) +  #NES:标准化ES，＜1.5下调
    geom_col(aes(fill=NES)) + #aes(fill= NES < 1.5) NES pval
    scale_fill_distiller(palette = "PiYG",direction = -1) + #RdPu 浅紫 YlOrRd黄色 PiYG红绿
    coord_flip() +
    labs(x="", y="Normalized Enrichment Score",
         title="Hallmark pathways") +
    theme_minimal() 
  fp
  #：https://blog.csdn.net/weixin_54199212/article/details/123536464
  ggsave("./out/GSEA_all_hallmark-barplot.pdf",width = 5,height = 4,family="serif")
  
  #5 HALLMARK-ssGSEA-group-pheatmap----
  
  geneset <- read.gmt("./data/h.all.v7.5.1.symbols.gmt")
  geneset <- split(geneset$gene, geneset$term)
  genesetInfo <- read.delim("./data/GenesetInfo.txt", sep = ",")
  genesetInfo$Classification <-ifelse(genesetInfo$Classification=="","Unclassified",genesetInfo$Classification)
  #genesetInfo <- subset(genesetInfo, Classification != "")
  levels = c("Immune", "Metabolism", "Signaling", "Proliferation","Unclassified")
  genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
  genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
  genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
  genesetInfo$Pathway <-gsub("HALLMARK_", "",genesetInfo$geneset )
  #genesetInfo$Pathway <-tolower(gsub("HALLMARK_", "", levels(genesetInfo$geneset)))
  library(GSVA)
  hallmark50_score <- gsva(as.matrix(expr[,-1]), #
                           geneset, #基因集
                           method = "ssgsea", kcdf = "Poisson", min.sz = 10)#
  #5.1差异热图#----
  library(limma)
  library(edgeR)
  #分组信息
  # group <-data.frame(name=c("C1","C2","C3","D1","D2","D3"),
  #                    group=rep(c('C', 'D'), each = 3, length.out = 6) )
  #group_list <-c("C1","C2","C3","D1","D2","D3")#rep(c('C', 'D'), each = 3, length.out = 6) #A vs B 组 谁是给药组？  B vs A  #c(rep("A",3),rep("B",3))
  #dgelist <- DGEList(counts = a549_exp, group = group)
  
  exprSet <- hallmark50_score #log2 ?
  #contra <-c('D-C')
  
  
  # 构建分组矩阵--design 
  design <- model.matrix(~0+factor(group$group))
  colnames(design) <- levels(factor(group$group))
  rownames(design) <- colnames(exprSet)
  contrast.matrix <- makeContrasts(contra,levels = design)
  fit <- lmFit(exprSet,design) #非线性最小二乘法
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
  DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
  DEG <- na.omit(DEG)
  DEG$Pathway <-rownames(DEG)
  # 输出差异分析表格
  write.table(DEG, file = "./out/output_hallmark.txt",
              sep = "\t", row.names = F, col.names = T, quote = F)
  
  #绘图#
  # 准备画图数据
  DPs <- read.table("./out/output_hallmark.txt", header = T)
  plot.data <- DPs
  plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
  plot.data$Celltypes <- contra #"IB vs Control"#factor(plot.data$Celltypes)
  plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset) # 通路顺序与genesetInfo一致
  plot.data$FDR<- cut(plot.data$adj.P.Val, breaks = c(0, 0.05, 0.5, 1),
                      include.lowest = T)
  plot.data$FDR <- factor(as.character(plot.data$FDR),
                          levels = rev(levels(plot.data$FDR)))
  levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))
  
  #
  library(ggplot2)
  color = c("#4682B4", "#FFFFFF", "#CD2626") # logFC对应的颜色
  #class.color = c("Immune" = "#D58986", "Metabolism" = "#80554C","Signaling" = "#71AC7A", "Proliferation" = "#E8D4B4") # 通路分类对应的颜色
  class.color = c("Immune" = "#A6CEE3", "Metabolism" = "#CAB2D6","Signaling" = "#B2DF8A", "Proliferation" = "#FDBF6F") # 通路分类对应的颜色
  #library(RColorBrewer);display.brewer.all();brewer.pal(12, "Paired")
  #"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00""#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
  
  p1 <- ggplot(plot.data, aes(x = Pathway, y = Celltypes, color = logFC, size = FDR#adj.P.Val#
  )) +
    geom_point() +
    scale_color_gradient2(low = color[1], mid = color[2], high = color[3]) +
    # geom_hline(yintercept = seq(min(as.numeric(plot.data$Celltypes))-0.5,
    #                             max(as.numeric(plot.data$Celltypes))+0.5),
    #            color = "grey80") +
    # geom_vline(xintercept = seq(min(as.numeric(plot.data$Pathway))-0.5,
    #                             max(as.numeric(plot.data$Pathway))+0.5),
    #            color = "grey80") +
    theme_classic() +
    theme(axis.title.y = element_blank(),axis.title = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.line = element_blank())+ #, legend.position = "top"
    coord_flip()
  p1
  ggsave("./out/hallmark_p1.pdf", width = 4, height = 8,family="serif")
  #ggsave(p1,filename="./out/hallmark_p1-1.png", width = 4, height = 8)
  
  p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) +
    geom_tile() +
    theme_classic() +
    scale_fill_manual(values = class.color) +
    theme(axis.text = element_blank(), axis.title = element_blank(),legend.position = "left",
          axis.ticks = element_blank(), axis.line = element_blank())+
    coord_flip()
  p2
  ggsave("./out/hallmark_p2_1.pdf", width = 1.55, height = 8,family="serif")
  
  library(cowplot)
  Sys.setenv(LANGUAGE = "en") #显示英文报错信息
  options(stringsAsFactors = FALSE) #
  #plot_grid(p1, p2, ncol = 1, align = 'v', rel_heights = c(10, 1)) # 拼接图像
  plot_grid(p2, p1, ncol = 2, align = 'h', rel_widths  = c(1, 3)) # 拼接图像
  ggsave("./out/hallmark_1.pdf", width = 6, height = 8,family="serif")
  
  #5.2 hallmark分组热图----
  # library(pheatmap);library(RColorBrewer)
  # pheatmap(hallmark50_score,
  #          #scale="column",#row
  #          cluster_rows = T,cluster_cols  = F, 
  #          fontsize_row = 8,fontsize_col = 10, #5 10
  #          show_colnames = T,show_rownames = T,
  #          border_color = "white",#cellwidth=10,cellheight = 10,
  #          gaps_col= 3, #table(exp_marker1$response)
  #          color= rev(colorRampPalette(brewer.pal(11, "PiYG"))(8))#Spectral,#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
  #          #annotation_names_row = F,annotation_row = anno_row_p,anno_names_row = F,annotation_colors = ann_colors_p
  # )
  # #display.brewer.all() PRGn 绿紫 PiYG酒红-浅绿
  #热图score
  score_df<-as.data.frame(hallmark50_score)
  rownames(score_df)<-gsub("HALLMARK_", "", rownames(score_df) )
  
  # score_df$Pathway<-rownames(score_df)
  # #score_df$Pathway <- factor(rownames(score_df), levels = genesetInfo$geneset) # 通路顺序与genesetInfo一致
  # #levels(score_df$Pathway) <- gsub("HALLMARK_", "", levels(score_df$Pathway))#tolower(gsub("HALLMARK_", "", levels(score_df$Pathway)))
  # score_df$Pathway <- gsub("HALLMARK_", "", score_df$Pathway)
  # #score_df<-na.omit(score_df)  #删除NA通路  45
  # 
  # score_df[,1:6] <-scale(score_df[,1:6])
  # rownames(score_df) <-score_df$Pathway
  # #排序
  # df_p <- t(score_df[,1:6])#t(scale( )) #Z-score
  # library (reshape2);library(dplyr)#数据转换
  # df_p <-df_p %>% melt() #长列表
  # colnames(df_p )<-c("sample","pathway","value")#sample pathway value
  # 
  # ggplot(df_p,aes(x=sample,y=pathway,fill=value,group=sample))+
  #   geom_tile(color="gray")+labs(x="",y="",fill="Score")+
  #   theme(legend.key.size = unit(0.15, "inches"),legend.title = element_text(size = 10),
  #         #axis.text.x =element_blank(),
  #         axis.ticks.x = element_blank(),panel.background = element_blank(),panel.grid = element_blank() )+ # + theme_bw()
  #   #scale_fill_gradientn(values = seq(0,1,0.2),colours = c('#6699CC','#FFFF99','#CC3333'))
  #   #scale_fill_viridis_c(option = "D",alpha=0.9)   #+facet_grid(~sample, scales = "free") #一模一样？组见差距较小
  #   scale_fill_distiller(palette = "Spectral") #scale_fill_gradientn(colors = brewer.pal(3, "Blues"))#scale_fill_brewer(palette="Set1")
  
  library(RColorBrewer)
  #display.brewer.all(type = "all")
  #display.brewer.all(type = "seq")
  
  library(pheatmap)
  pheatmap(score_df[,1:6])
  #rownames(score_df)<-gsub("HALLMARK_", "",rownames(score_df))
  
  annotation_col <-data.frame(
    Group=group$group )
  rownames(annotation_col) <-group$name#colnames(ferro_score1[,1:6])
  
  ann_colors <-list(Type=c("Immune" = "#A6CEE3", "Metabolism" = "#CAB2D6","Signaling" = "#B2DF8A", "Proliferation" = "#FDBF6F","Unclassified"="grey70") #,  #'Unchange'='#CAB2D6',
                    #Group=c(paste0(unique(group$group)[1]) ="#7EC7A7",paste(unique(group$group)[2]) ="skyblue") 
  )
  
  anno_row <- data.frame(
    Type=genesetInfo$Classification) #, size = 20, replace = TRUE  adj.P.Val
  rownames(anno_row) <- gsub("HALLMARK_", "", genesetInfo$Pathway)#gsub("HALLMARK_", "", levels(genesetInfo$geneset)) #tolower(gsub("HALLMARK_", "", levels(genesetInfo$geneset)))
  
  #plot
  hallmark_p <-pheatmap(score_df[,1:6],cluster_rows = T,cluster_cols  = F, show_colnames = T,
                        annotation_names_row = F, annotation_names_col = F,
                        border_color = "white",#cellwidth=10,cellheight = 10,
                        gaps_col= 3, #table(exp_marker1$response)
                        color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)) ,#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
                        annotation_col=annotation_col,annotation_colors = ann_colors,
                        annotation_row = anno_row
  )
  hallmark_p
  save_pheatmap_pdf <- function(x, filename, width, height) {
    library(grid)
    x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.02)#修改聚类树线条宽度
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height,family = "serif")
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }#保存热图的函数
  save_pheatmap_pdf(hallmark_p, "./out/hallmark_pheatmap-1.pdf",7,7) #6 6
  
  ##添加显著性
  pvalue_g<-data.frame() #T.test
  #score_df<-as.data.frame(hallmark50_score)
  for (i in 1:nrow(score_df) ){
    print(i)
    pwilcox <- t.test(as.numeric(score_df[i,4:6] ),as.numeric(score_df[i,1:3]) )
    fc <- mean(as.numeric(score_df[i,4:6] ) )/mean( as.numeric(score_df[i,1:3]) )#drug/control
    pvalue_g<-rbind(pvalue_g,data.frame(pvalue=pwilcox$p.value,FoldChange=fc,gene= rownames(score_df)[i]) )
  } #wilcox.test 均不显著；t.test 较多显著性!
  #对样本scale后计算,0则NA
  
  pvalue_g$Sign. <- ifelse(pvalue_g$pvalue>= 0.05,"ns",
                           ifelse(pvalue_g$pvalue < 0.01,"**","*") )
  pvalue_g$FC <-ifelse(pvalue_g$FoldChange <=1,"< 1",
                       ifelse(pvalue_g$FoldChange >2,"> 2","> 1") )
  rownames(pvalue_g) <-pvalue_g$gene
  
  annotation_row_p <- data.frame(
    Sign.=pvalue_g$Sign.,FoldChange=pvalue_g$FC) #, size = 20, replace = TRUE  adj.P.Val
  rownames(annotation_row_p) <- rownames(pvalue_g)
  ann_colors_p=list(FoldChange=c('< 1'='#CAB2D6','> 1'="#FDBF6F",'> 2'="#B15928"),
                    Sign.=c("ns"="grey","*"="#B2DF8A",'**'="#33A02C" ))
  
  library("gplots");library(viridis)
  hallmark_g <- pheatmap(score_df, #exp_genes[rownames(exp_genes) !=c("CD44","CD9"),],
                         scale="column",border_color = "white",#NA,
                         cluster_rows = T,cluster_cols  = F,gaps_col= 3,
                         annotation_row = annotation_row_p,annotation_names_row = F,
                         annotation_colors = ann_colors_p,annotation_col=annotation_col,
                         annotation_names_col = F,
                         color= rev(colorRampPalette(brewer.pal(10, "Spectral"))(15) ),
                         #color=viridis(20, option = "G")[3:20], #G D
                         #color=bluered(50)#library("gplots")
  )
  hallmark_g
  save_pheatmap_pdf(hallmark_g, "./out/hallmark_pheatmap-p-fc.pdf",7,7) 
  
}

enrichment(gene_type,input_gene,exprSet,group,contra)


#3.0 Drug vs. Control T.test----
normal <-group[group$group=="NC",]$name #colnames(exp)[grep("P", colnames(exp) )] #paired normal 
drug <- group[group$group=="OE",]$name#colnames(exp)[grep("T", colnames(exp) )] #tumor 

pvalue_g<-data.frame()
for (i in rownames(exp) ){ #[1:10]
  print(i)
  pwilcox <- t.test(as.numeric(exp[i,drug]),  #T
                         as.numeric(exp[i,normal]) ) #P
  fc <- mean( as.numeric(exp[i,drug]) )/mean( as.numeric(exp[i,normal]) )#drug/control
  #fc <- as.numeric(pwilcox$estimate[1])/as.numeric(pwilcox$estimate[2])
  pvalue_g<-rbind(pvalue_g,data.frame(pvalue=pwilcox$p.value,FoldChange=fc,gene=i  ) )
  
}
pvalue_g<-na.omit(pvalue_g)
pvalue_g$Sign. <- ifelse(pvalue_g$pvalue>= 0.05,"ns",
                         ifelse(pvalue_g$pvalue < 0.01,"**","*") )
pvalue_g$Fold <-ifelse(pvalue_g$FoldChange <=1,"< 1",
                       ifelse(pvalue_g$FoldChange >2,"> 2","> 1") )
rownames(pvalue_g) <-pvalue_g$gene
write.csv(pvalue_g,"./out/OE-STAT3-OE-NC_ttest.csv")