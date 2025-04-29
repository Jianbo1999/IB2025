#ROI tumor only


if(T){
  rm(list = ls())
  setwd("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\R\\1.DEG_enrichment\\1.1_Tumor")
  
  #输出文件夹!!!!!!!
  program_name="IB_DM_Tumor"
  
  folder_path <- program_name#"./out" # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
}
load("G:/IB/231-202306/4T1动物/DSP/测序结果/R/1.DEG_enrichment/1.1_Tumor/1.1_IB_DM_Tumor.RData")

#0. load exp and meta----
if(T){
  exp<-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\1.Quality_Control\\raw_exp.csv",row.names = 1)
  meta<-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\1.Quality_Control\\overall_annotation.csv",row.names = 1)
  meta$id <-rownames(meta)#ID
  #log2  表达矩阵
  exprSet<-log2(exp + 1)   #
  # boxplot(exp[1:10])
  # boxplot(exprSet[1:10])#5-15
  ####tumor only
  meta<-meta[meta$Type=="TC",]#30
  exprSet<-exprSet[,colnames(exprSet) %in% meta$id ] #30 ROI
  #
  meta$group <-meta$Group #比较分组信息!!!!!!!!!!
  group <-meta[,c("id","group")]
  contra <-c('IB-DM')#对比信息
  adjPval <- 0.01         # 矫正p值-p阈值
  aflogFC <- 1            # logFC  阈值
}
#0.1 pca----
#LOG2 +1
library(ropls);library(FactoMineR) # PCA函数
library(factoextra) # fviz_pca_ind函数

###PCA分析 
data.pca <- prcomp(t(exprSet), scale. = T)  #这是后续作图的文件 
#输出特征向量 
write.csv(data.pca$rotation, file=paste0(program_name,"/PCA_IB_DM.csv") ) #输出新表 write.table(predict(data.pca), file="newTab.xls", quote=F, sep = "\t") #输出PC比重 pca.sum=summary(data.pca) write.table(pca.sum$importance, file="importance.xls", quote=F, sep = "\t") 

##plot
library(factoextra)
fviz_pca_ind(data.pca, 
             col.ind=group$group,  ####drug
             #col.ind=meta$Type,####PCA-location
             geom = c("point"),
             mean.point=F, addEllipses = F, #ellipse.level=0.95,ellipse.type="confidence", 
             title="",legend.title="Group", 
             palette = c("skyblue", "pink"))+  #"skyblue", "pink","#7EC7A7"
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.1, linetype="solid"))#加个边框
ggsave(paste0(program_name,"/PCA_IB_DM_0.pdf"),width = 3.5,height = 3,family="serif")

fviz_pca_ind(data.pca, 
             col.ind=group$group,  ####drug
             #col.ind=meta$Type,####PCA-location
             geom = c("point"),
             mean.point=F, addEllipses = T, #ellipse.level=0.95,ellipse.type="confidence", 
             title="",legend.title="Group", 
             palette = c("skyblue", "pink"))+  #"skyblue", "pink","#7EC7A7"
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.1, linetype="solid"))#加个边框

ggsave(paste0(program_name,"/PCA_IB_DM.pdf"),width = 3.5,height = 3,family="serif")


#勉强
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
  #筛选差异表达基因 adjPval差
  DEG$Group <- ifelse((DEG$P.Value <adjPval & abs(DEG$logFC)>aflogFC), ifelse(DEG$logFC>aflogFC,"Up","Down"), "Stable")
  print(table(DEG$Group))
  # 输出差异分析表格
  write.csv(DEG,paste0(program_name,"/all_DEG_limma.csv") )#down 272 stable 19736
  #write.table(DEG, file = "./program_name/all_DEG_limma.txt",sep = "\t", row.names = F, col.names = T, quote = F)
  
}
# Down Stable 
# 305  19703 
#重新设置阈值？0.5 --> -1？


##volcano----
library(dplyr)
top_5 <-DEG %>% arrange(P.Value,desc(abs(logFC)) ) %>% head(5)
gene_volcano<-top_5$Symbol
#gene_volcano <-  DEG[abs(DEG$logFC) >1,]$Symbol #DEG[DEG$Group != "Stable",]$Symbol#DEG[DEG$Group != "Stable",]$Symbol#rownames(DEG[abs(DEG$logFC) >1.5,])#火山图展示基因 
gene_volcano <-c(gene_volcano,"Cd274") #输入兴趣基因

#c("STAT3","STAT4","SLC47A1","SLC7A11","GPX4","CD274")#

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

p1 #丑
ggsave(paste0(program_name,"/volcano-D-1.pdf"),family="serif",width = 5,height = 4)

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
                  size = 3, 
                  #box.padding = 0.5, #字到点的距离 
                  #point.padding = 0.8, #字到点的距离，点周围的空白宽度 
                  min.segment.length = 0.1, 
                  segment.size=0.2, 
                  color = 'black') +
  labs(x="Log2(FC)",y="-Log10(pvalue)") #adj.Pval
ggsave(paste0(program_name,"/volcano-2.pdf"),family="serif",width = 5,height = 4)

##ggVolcano
#https://github.com/BioSenior/ggVolcano

library(ggVolcano)
#deg_data<-deg_data
log2FC_cut=1;FDR_cut=0.01

data <- add_regulate(DEG, log2FC_name = "logFC",
                     fdr_name = "P.Value",log2FC = log2FC_cut, fdr = FDR_cut) #无上调基因无法使用？
data$row <-data$Symbol

ggvolcano(data, x = "log2FoldChange", y = "padj",legend_position="UR",
          log2FC_cut=log2FC_cut,FDR_cut=FDR_cut,custom_label=gene_volcano, #兴趣基因
          x_lab="Log2 (fold change)",y_lab="-Log10 (p value)",
          fills = c("skyblue","grey90","pink"),
          colors = c("skyblue","grey90","pink"),
          label = "row", label_number = 10, output = FALSE)
ggsave(paste0(program_name,"/ggvolcano.pdf"),family="serif",width = 4,height = 4)
## gradient color
gradual_volcano(data, 
                log2FC_cut=log2FC_cut,FDR_cut=FDR_cut,custom_label=gene_volcano,
                x_lab="Log2 (fold change)",y_lab="-Log10 (p value)",legend_title="-Log10 (p value)",
                x = "log2FoldChange", y = "padj",legend_position="UR",
                label = "row", label_number = 10, output = FALSE)+
  ggsci::scale_color_gsea()+ggsci::scale_fill_gsea()
ggsave(paste0(program_name,"/gradual_volcano.pdf"),family="serif",width = 4,height = 4)

##GO term volcano plot
# library(RColorBrewer)
# # Change the fill and color manually:
# deg_point_fill <- brewer.pal(5, "RdYlBu")
# names(deg_point_fill) <- unique(term_data$term)
# data("term_data") #需要输入富集信息！！Gene.names term
# term_volcano(data, term_data,
#              x = "log2FoldChange", y = "padj",
#              normal_point_color = "#75aadb",
#              deg_point_fill = deg_point_fill,
#              deg_point_color = "grey",
#              legend_background_fill = "#deeffc",
#              label = "row", label_number = 10, output = FALSE)
##rank----
#https://mp.weixin.qq.com/s/Go1_K6MtONJyhkpbmm69Bg

if(T){
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(ggfun)
  
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
  
  ggsave(filename = paste0(program_name,"/DEG_rank.pdf"),
         height = 7,
         width = 9)
}

##兴趣基因热图heatmap----

library(pheatmap);library(RColorBrewer)
#pheatmap(exprSet)#太多基因，卡
#差异基因 DEG[DEG$Group != "Stable",]$Symbol
#gene_volcano1 =gene_volcano
gene_volcano1 =DEG[DEG$Group != "Stable",]$Symbol #272
hallmark_p <-pheatmap(exprSet[rownames(exprSet) %in% gene_volcano1, ], #gene_volcano
                      #exprSet[rownames(exprSet) %in% c(DEG[DEG$Group != "Stable",]$Symbol) , ]##差异基因 
                      scale="column",#row
                      cluster_rows = T,cluster_cols  = F, 
                      fontsize_row = 8,fontsize_col = 10, #5 10
                      show_colnames = F,show_rownames = T,
                      border_color = "white",#cellwidth=10,cellheight = 10,
                      #gaps_col= 3, #table(exp_marker1$response)
                      color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)) ,#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
                      annotation_col=meta[c(32,31,33)]#,annotation_row = anno_row
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
save_pheatmap_pdf(hallmark_p, paste0(program_name,"/DEG_genes_pheatmap.pdf"),6,18) 

##某个基因T检验boxplot----
gene <- "Cd274"

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
    geom_point(alpha=0.5,size=1.5,
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
    stat_compare_means(method = "wilcox.test",size=3.5,show.legend= F,label = "p.format",
                       label.x =1.2,label.y =max(df_gene$Gene)-max(df_gene$Gene)/20 )#+ #,label.y =max(df_gene$Gene)-0.55,label.x =1.5,label.y = 1.3  0.8 ,label.y = 4 #非参数检验kruskal.test not anova
  #ylim(NA,max(df_gene$Gene)+0.5) #wilcox.test  t.test
  #ggsave(paste0("./out/boxplot_Ttest_",gene,".pdf"),width = 2,height = 3,family="serif")
  ggsave(paste0(program_name,"/boxplot_",gene,"_wilcox.pdf"),width = 2,height = 3,family="serif")
  
}

##批量基因boxplot----
gene_more <-gene_volcano1 #unique(c(gene_volcano,c("STAT3","STAT4","SLC47A1","SLC7A11","GPX4","CD274"))) #输入基因集

for (i in gene_more ){
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
    geom_point(alpha=0.5,size=1.5,
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
    stat_compare_means(method = "wilcox.test",size=3.5,show.legend= F,label = "p.format",
                       label.x =1.2,label.y =max(df_gene$Gene)-max(df_gene$Gene)/20 )
  ggsave(paste0(program_name,"/boxplot_",i,"_wilcox.pdf"),width = 2,height = 3,family="serif")
  
}

##gene_volcano top 热图+显著性----
#T检验显著性表格
if(T){
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
  col<-data.frame(group$group);rownames(col)=rownames(group);colnames(col)="Group"
  int_p <-pheatmap(exp_int,
                   scale="column",#row
                   cluster_rows = T,cluster_cols  = F, 
                   fontsize_row = 8,fontsize_col = 10, #5 10
                   show_colnames = F,show_rownames = T,
                   border_color = "white",#cellwidth=10,cellheight = 10,
                   #gaps_col= 3, #table(exp_marker1$response)
                   color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)),#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
                   annotation_col=col,annotation_names_row = F,annotation_row = anno_row_p,anno_names_row = F,annotation_colors = ann_colors_p
  )
  int_p
  save_pheatmap_pdf(int_p, paste0(program_name,"/int_pheatmap-1.pdf"),5,4) #allDEGgene
  
}

#2.0 富集分析----
if(T){
  ########
  #使用差异基因
  input_gene <-DEG[DEG$Group != "Stable",]$Symbol #272
  gene_type<-"ALL" #上调基因？
  OrgDb = 'org.Mm.eg.db'#小鼠 #"org.Hs.eg.db" #人
  KEGG_database <- 'mmu'
  #data文件夹paste0(program_name,"/int_pheatmap-1.pdf")
  
  #2.1.0 GO----
  library(clusterProfiler)
  trans = bitr(input_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  library(dplyr)
  #gene_df <- data.frame(gene = trans$ENTREZID)
  # eego <- enrichGO(gene          = trans$ENTREZID,
  #                  #universe      = names(geneList),
  #                  OrgDb         = OrgDb,keyType = "ENSEMBL",
  #                  ont           = "ALL",
  #                  pAdjustMethod = "BH",
  #                  #pvalueCutoff  = 0.01,
  #                  qvalueCutoff  = 0.05,
  #                  readable      = TRUE)
  eego<-enrichGO( trans$ENTREZID,#GO富集分析
                  OrgDb = OrgDb,
                  keyType = "ENTREZID",#设定读取的gene ID类型
                  ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                  pvalueCutoff = 1,#设定p值阈值
                  qvalueCutoff = 1,#设定q值阈值
                  readable = T)
  #plot
  # df_go <- data.frame(eego) %>%
  #   group_by(ONTOLOGY) %>%
  #   slice_head(n = 10) %>% #前10个,前5个？
  #   arrange(desc(pvalue))
  df_go<-eego@result%>% #slice_head(n = 10) %>% #前10个,前5个？
    arrange(desc(Count)) #-desc(pvalue)
  ratio <- lapply(df_go$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
  df_go$ratio <- ratio
  df_go$Description <- factor(df_go$Description,levels = df_go$Description)
  df_go<-df_go[df_go$pvalue <0.05,]
  write.csv(df_go,paste0(program_name,"/DEG_",gene_type,"_GO.csv"))
  
  #plot
  my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  library(ggplot2)
  ## 选择部分通路展示！！！！！！！！！
  #dotplot(eego,showCategory = 5)#+facet_grid(ONTOLOGY~., scale="free")#点状图
  dotplot(eego,color = "pvalue", #p.adjust title = 'GO',
          showCategory = c(df_go$Description[c(8,11,32,72,253,325,330,390,418,479,552,561,747,783,908)] )
  )+
    scale_color_continuous(low='pink', high='#57C3F3')
  ggsave(paste0(program_name,"/DEG_",gene_type,"_GO_my-1.pdf"),width = 5.5,height = 7,family="serif")
  barplot(eego,showCategory = 10,title = 'GO')
  ggsave(paste0(program_name,"/DEG_",gene_type,"_GO_barplot10.pdf"),width = 5.5,height = 6,family="serif")
  
  ##GGPLOT2-barplot
  #df$percentage <- paste0(round(df$value / sum(df$value) * 100, 1), "%")

  for (i in c("RdPu","YlOrRd","PiYG") ){
    print(i)
    ggplot(data = df_go[c(8,11,32,72,253,325,330,390,418,479,552,561,747,783,908),], 
           aes(x = Count, y = Description, fill =-log10(pvalue) )) +  #Description, fill = Cluster
      #scale_fill_viridis_c(alpha = 0.7,option = "G")+ #,direction = -1
      #scale_fill_manual(values =rev( colorRampPalette(brewer.pal(11, "Paired"))(length(kkgsea_res_p$Description)) )  ) +
      scale_fill_distiller(palette = i,direction = 1) + #RdPu YlOrRd PiYG
      geom_bar(stat = "identity", width = 0.7, alpha = 0.5) +
      #scale_x_continuous(expand = c(0,0),limits = c(min(kkgsea_res_p$NES)-0.5,max(kkgsea_res_p$NES)+0.5)) + # 调整柱子底部与y轴紧贴
      #scale_y_continuous(expand = c(0,-1)) +
      labs(x = "Count", y = "", title = "") +
      geom_text(aes(x=0.03, # 增大与条形图的间距
                    label = Description,
                    hjust = 0 ),
                size = 3.5) +
      theme_classic() +theme(axis.text.y = element_blank(),axis.ticks.y = element_blank() )
    ggsave(paste0(program_name,"/DEG_",gene_type,"_GO_ggplot2_",i,".pdf"),family="serif",width = 4.5,height = 4)
    
  }
  
  #2.2.0 KEGG----
  kegg <- enrichKEGG(
    gene          = trans$ENTREZID,
    keyType     = "kegg",
    organism   = KEGG_database,#'hsa',
    pvalueCutoff      = 0.05,
    pAdjustMethod     = "BH",
    qvalueCutoff  = 0.05
  )
  library(DOSE)
  kegg<-setReadable(kegg, OrgDb = 'org.Mm.eg.db',#"org.Hs.eg.db", 
                    keyType="ENTREZID")
  
  #2.1 enrichplot
  logFC <- DEG[ DEG$Symbol %in% input_gene,]$logFC
  names(logFC) <- DEG[ DEG$Symbol %in% input_gene,]$Symbol  
  
  library(enrichplot)
  
  pdf(paste0(program_name,"/DEG_",gene_type,"_KEGG_cnetplot.pdf"),family="serif",width = 12,height = 4)
  enrichplot::cnetplot(kegg,showCategory = 10,#foldChange = logFC,
                       circular=F,colorEdge = TRUE)
  enrichplot::cnetplot(kegg,showCategory = 10,foldChange = logFC,
                       circular=F,colorEdge = TRUE)
  # cnetplot(kegg, showCategory = 10, #选择top10的pathway ，这里也可以用包含pathway名称的向量           
  #          color.params = list(foldChange = logFC, #用上面的logFC值标注基因的颜色
  #                                   edge = T)) #显示pathway的标注
  # heatplot(kegg, foldChange=logFC,
  #          #showCategory = 8,
  #          symbol = "dot", #"rect" 矩形
  #          pvalue = NULL,
  #          label_format = 30)
  dev.off()
  
  ego <- kegg@result
  ego <- ego [order(ego$Count,decreasing = T),]
  ego <- ego[ego$pvalue <=0.05,] #P > 0.05
  #ego <- ego[,c(2,5,8,9)]#
  
  library(ggplot2);library(ggsci);library("scales")
  ego$Description <- factor(ego$Description,levels = c(ego$Description))
  ego$Description<-gsub(" -.*","",ego$Description )
  write.csv(ego,paste0(program_name,"/DEG_",gene_type,"_KEGG.csv"))
  #选择通路！！！！！！！！！！！！兴趣通路-见CSV
  ego1<-ego[c(22,37,46,50,57,66,70,73,81,109,152),]
  
  #2.2 barplot
  mytheme<- theme(axis.title = element_text(size = 13),
                  axis.text = element_text(size = 11),
                  plot.title = element_text(size = 14,
                                            hjust= 0.5,
                                            face= "bold"),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 11))
  #批量颜色
  for (i in c("RdPu","YlOrRd","PiYG") ){
    print(i)
    ggplot(data = ego1,
           aes(x = Count, y = reorder(Description,-Count), fill =-log10(pvalue) )) +  #Description, fill = Cluster
      scale_fill_distiller(palette = i,direction = -1) + #RdPu YlOrRd PiYG#scale_fill_distiller(palette = i,direction = 1) + #RdPu YlOrRd PiYG
      geom_bar(stat = "identity", width = 0.7, alpha = 0.5) +
      labs(x = "Count", y = "", title = "KEGG") +
      geom_text(aes(x=0.03, # 增大与条形图的间距
                    label = Description,
                    hjust = 0 ),
                size = 3.5) +
      theme_classic() +theme(axis.text.y = element_blank(),axis.ticks.y = element_blank() )
    ggsave(paste0(program_name,"/DEG_",gene_type,"_KEGG_bar_",i,".pdf"),family="serif",width = 5,height = 4)
    
    #2.3 dotplot
    p<- ggplot(data = ego1,
               aes(x = Count, y = reorder(Description,-Count) ) )+ #, fill = -log10(pvalue)
      scale_color_distiller(palette = i,direction = -1) + #RdPu YlOrRd PiYG
      ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 
      #geom_bar(stat = "identity", width = 0.8) +
      geom_point(aes(size=Count,color=-log10(pvalue) ))+
      labs(x = "Count", y = "", title = "KEGG") + 
      theme_bw()+ mytheme+
      scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) #Y axis next line换行
    #scale_color_viridis_c(alpha = 0.7)
    p
    ggsave(paste0(program_name,"/DEG_",gene_type,"_KEGG_dotplot_",i,".pdf"),family="serif",width = 5.5,height = 3.5)
    
      
  }
  
 
  #2.4 barplot+基因
  #https://mp.weixin.qq.com/s/vZXxoXSOX6Y7bnIuTHaBzg
  # 先自定义主题： #若基因太大则无法显示完全
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
  p <- ggplot(data = ego1, aes(x = -log10(pvalue), y = reorder(Description,-Count), fill =Description )) +  #, fill = Cluster
    #scale_fill_manual(values =c('#6bb9d2', '#d55640')) +
    scale_fill_manual(values =my36colors) +
    #scale_fill_manual(values =rev( colorRampPalette(brewer.pal(11, "Paired"))(length(ego$Description)) )  ) +
    geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
    scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
    #scale_y_continuous(expand = c(0,-1)) +
    labs(x = "-log10(p value)", y = "Pathway", title = "KEGG") +
    # x = 0.61 用数值向量控制文本标签起始位置
    geom_text(size=3.8, aes(x = 0.05, label = Description),hjust = 0) + # hjust = 0,左对齐
    #geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust =2.5, color="skyblue" ) + # hjust = 0,左对齐
    geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust =-1, color="skyblue" ) + # hjust = 0,左对齐
    theme_classic() + mytheme1 
  #ylim(min(ego$Description),max(ego$Description))
  #scale_y_discrete(labels = \(x) sub(pattern = "_", replacement = " ", x, fixed = TRUE)) 
  p
  ggsave(paste0(program_name,"/DEG_",gene_type,"_KEGG_bar+genes.pdf"),width = 6,height = 3.5,family="serif")#7.5
  
  #2.3.0 GSEA----
  
  #输入所有基因！！！
  trans_all = bitr(DEG$Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb) #"org.Hs.eg.db"
  colnames(trans_all)[1]<-"Symbol"
  allDEG <-merge(DEG,trans_all,by="Symbol")
  
  genelist <-  allDEG$logFC#所有基因-18866
  names(genelist) = allDEG$ENTREZID
  genelist = sort(genelist, decreasing = TRUE)#排序
  any(genelist == 0)
  genelist<-genelist[genelist>0] #Error in if (abs(max.ES) > abs(min.ES)) { :    missing value where TRUE/FALSE needed
  kkgsea <- gseKEGG(geneList     = genelist,
                    organism     = KEGG_database,#'hsa', 
                    minGSSize    = 1,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.9,#use_internal_data=T,
                    pAdjustMethod = "none" ) #进行gseKEGG富集分析 
  kkgsea=setReadable(kkgsea,keyType = 'ENTREZID',OrgDb = OrgDb) #'org.Hs.eg.db'
  kkgsea_res <- kkgsea@result 
  kkgsea_res_p <- kkgsea_res[kkgsea_res$p.adjust <0.05,]#
  write.csv(kkgsea_res, paste0(program_name,"/GSEA_all_gene.csv") ) #paste0("./out/gsea_",gene_type,".csv")
  #无显著性结果！
  
  library(enrichplot)
  #gseaplot2(kkgsea,geneSetID ="hsa00310",color = "red" ) #
  pdf(paste0(program_name,"/GSEA_enrichplot_top10.pdf"),width = 6,height = 3.5,family="serif")
  gseaplot2(
    kkgsea, #gseaResult object，即GSEA结果
    #geneSetID = c(1:10),#"hsa00310",#1,#"hsa05168",#富集的ID编号
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
  #ggsave("./out/GSEA_all_barplot.pdf",width = 6.5,height = 4)
  
  ##3.2 GSEAvis----
  #3.2.1 gseaNb
  library(GseaVis)
  #批量输出所有显著性KEGG
  pdf(paste0(program_name,"/GSEA_all_GSEAVIS.pdf"),width = 4,height = 4,family = "serif")
  for (i in kkgsea_res$ID){ #kkgsea_res_p
    print(i)
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
  
  
  # #expr <-readRDS("./out/exp_all_fpkm.rdata")[,7:12]
  # expr$gene_name <-rownames(expr)
  # expr<- expr %>% dplyr::select('gene_name',colnames(expr)[1:dim(expr)[2]-1],everything())
  # ##表达量文件第一列必须是基因名
  # 
  # 
  # pdf(paste0("./out/GSEA_all_GSEAVIS_pheatmap.pdf"),width = 6.5,height = 5.2)
  # for (i in kkgsea_res_p$ID){
  #   gsea_name <-as.character(
  #     kkgsea_res_p[kkgsea_res_p$ID ==i,]$Description
  #   )
  #   
  #   gsea_plot<-gseaNb(object = kkgsea,newGsea = T,
  #                     geneSetID = i, #gse@result[["ID"]][1]
  #                     add.geneExpHt = T,
  #                     exp = expr,
  #                     exp.col = c('skyblue','white','pink'),#修改热图颜色
  #                     ght.geneText.size = 10, #调整热图基因标签大小
  #                     ght.relHight = 0.3, #修改热图相对高度
  #                     addPval = T,#Add NES and Pvalue
  #                     pvalX = 0.6,pvalY = 0.7, #调整标签位置和颜色
  #                     pCol = 'red',pHjust = 0,
  #                     subPlot = 2  #只保留上2部分
  #   )
  #   print(gsea_plot)
  #   #ggsave(paste0("./out/GSEA_",gsea_name,".pdf"),width = 4,height = 4)
  # }
  # dev.off()
  
  #2.4.0 HALLMARK-GSEA----
  library(msigdbr) 
  hallmark_genesets <- msigdbr(species = "Mus musculus",#"Homo sapiens", 
                               category = "H")#下载基因集
  #table(hallmark_genesets$gs_name)
  
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
              filter(padj < 0.05) ,#%>%    #过滤筛选 pval
            #head(n= 20), 
            aes(reorder(pathway, NES), NES)) +  #NES:标准化ES，＜1.5下调
    geom_col(aes(fill=padj)) + #aes(fill= NES < 1.5) NES pval
    scale_fill_distiller(palette = "RdPu",direction = 1) + #RdPu 浅紫 YlOrRd黄色 PiYG红绿
    coord_flip() +
    labs(x="", y="Normalized Enrichment Score",
         title="Hallmark pathways") +
    theme_classic()#theme_minimal() 
  fp
  ggsave(paste0(program_name,"/GSEA_all_hallmark-barplot.pdf"),width = 5,height = 6,family="serif")
  
  
  ##HALLMARK GSEA图 
  #使用clusterProfiler中的read.gmt函数读取下载的gmt基因集文件
  
  gsea.hallmark<- GSEA(genelist_1,  #待富集的基因列表
                       TERM2GENE = hallmark_genesets[3:4],#fgsea_sets,  #基因集
                       pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                       pAdjustMethod = 'BH')  #指定 p 值校正方法
  gsea_res <- gsea.hallmark@result 
  gsea_res_p <- gsea_res[gsea_res$p.adjust < 0.05,]#
  write.csv(gsea_res, paste0(program_name,"/GSEA_all_gene_hallmark.csv") ) #paste0("./out/gsea_",gene_type,".csv")
  
  library(GseaVis) #批量输出所有显著性KEGG
  #pdf(paste0(program_name,"/GSEA_all_GSEAVIS_hallmark.pdf"),width = 4,height = 4,family = "serif")
  for (i in gsea_res_p$ID){ #kkgsea_res_p
    print(i)
    gsea_name <-as.character(
      gsea_res_p[gsea_res_p$ID ==i,]$Description
    )
    
    gsea_plot<-gseaNb(object = gsea.hallmark,
                      geneSetID = i, termWidth = 35,#hsa00310 hsa00480
                      #newGsea = F,addPoint = F,
                      #rmSegment = T, #移除红色标记线
                      #rmHt = T,#移除热图
                      addPval = T,#Add NES and Pvalue
                      pvalX = 0.35,pvalY = 0.6, #调整标签位置和颜色
                      pCol = 'red',pHjust = 0,subPlot = 2)
    print(gsea_plot)
    ggsave(paste0(program_name,"/GSEA_",gsea_name,".pdf"),width = 4,height = 4)
  }
  #dev.off()
  
  
  #2.6 通路活性评分progeny----
  library(progeny)
  pathways <- progeny(as.matrix(exprSet),  #matrix
                      scale=TRUE,
                      organism="Mouse", #如为小鼠，填 "Human" "Mouse"
                      top = 100,  #The top n genes for generating the model matrix according to significance (p-value)
                      perm = 1)
  pathways <- pathways %>% as.data.frame()%>% mutate(id=rownames(pathways)) #合并到mata
  pathways <-merge(pathways,meta,by="id") #2-15 14 pathways score
  library(tidyverse)
  pathways <-column_to_rownames(pathways,"id")
  
  ##pheatmap
  x <-select(pathways,group);rownames(x)<-rownames(pathways)
  ann_colors = list(
    group = c(DM="skyblue", IB="pink"),
    Type= c(IM=my36colors[1],TC=my36colors[2],TC_IM=my36colors[3]) )
  pathways <-arrange(pathways,group)
  
  library(pheatmap);library(RColorBrewer)
  pheat <-pheatmap(as.data.frame(t(pathways[,1:14])),show_colnames = F, #col sample; row genes
                   cluster_rows = T,cluster_cols = F,border_color = NA,gaps_col = as.numeric(table(pathways$group)[1]) ,
                   annotation_col=pathways[c(45,50)],#x, 
                   annotation_colors = ann_colors, 
                   legend=TRUE ,annotation_names_col=F,#legend_breaks=c(-2,0,2),  #设置图例范围
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlGnBu")))(10) ) # RdYlBu
  pheat
  #IB DM 效果差
  
  save_pheatmap_pdf(pheat, paste0(program_name,"/progeny_pheatmap_1.pdf"),4.5,3) 
  
  ##boxplot
  pathways_long <- pathways %>% tidyr::pivot_longer(cols = colnames(pathways)[1:14],
                                                    names_to = "pathway")
  library(ggplot2);library(ggsci);library(ggpubr)
  ggplot(data = pathways_long,
         aes(x = group, y = value, fill = group))+
    facet_wrap(pathway~. , scales='free',ncol = 7) +
    scale_fill_manual(values = c("#008ECB", "#EA921D", "#D14039")) + 
    geom_violin(alpha=0.4, position = position_dodge(width = .75),
                size=0.1, color="black") + # 边框线黑色
    geom_boxplot(notch = F, outlier.size = -1,
                 width=0.5,
                 color="black", lwd=0.5, alpha = 0.7)+ # 背景色透明化
    # geom_point(shape = 21, size=2, 
    #            position = position_jitterdodge(), 
    #            color="black", alpha=1)+ # 边框线黑色
    ylab(expression("score")) +
    xlab("")  +
    #ggtitle("TCGA")+
    stat_compare_means(method = "wilcox.test",
                       label.y = max(pathways_long$value),label.x = 1.5,
                       comparisons=list(c("IB","DM")),#comparisons=list(c("ClusterA","ClusterB"),c("ClusterB","ClusterC"),c("ClusterA","ClusterC")),
                       step.increase = 0.05,
                       label="p.signif",
                       vjust = 0.7
    )+#mytheme1+
    theme(legend.position='none')+theme_bw()
  ggsave(paste0(program_name,"/progeny_boxplot.pdf"),width = 12,height=6,family="serif")
  
  ##批量BOXPLOT
  for (i in colnames(pathways[1:14]) ){
    print(i)
    # df_gene <-exprSet[rownames(exprSet) %in% i,]
    # df_gene <-as.data.frame(t(df_gene))
    # df_gene <-as.data.frame(scale(df_gene) )  #scale样本
    # df_gene$Group <- group$group #c(rep("Control",3),rep("Drug",3))  #分组信息！！！
    
    df_gene<-pathways[,c(i,"Group")]
    colnames(df_gene)[1] <-"Gene"
    library(ggplot2);library(ggpubr)
    ggplot(df_gene,#
           aes(Group,Gene,color=Group))+  #age_group
      geom_boxplot(alpha=1,width=0.45,fill=NA,
                   position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      geom_point(alpha=0.5,size=1.5,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.8))+
      geom_violin(alpha=0.2,width=0.9,
                  position=position_dodge(width=0.8),
                  size=0.25)+
      scale_color_manual(values =  c("#EDA065","#66CCFF","#7EC7A7"))+ #blue c("#0066CC","#66CCFF","#9DC3E6")
      theme_bw() +labs(color="",title = i,y="Activity score",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
            panel.grid = element_blank(),legend.key = element_blank(),legend.position="none" ) +
      # stat_compare_means(#aes(group = Group) ,
      # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
      # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
      # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
      # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
      stat_compare_means(method = "wilcox.test",size=3.5,show.legend= F,label = "p.format",
                         label.x =1.2,label.y =max(df_gene$Gene)-max(df_gene$Gene)/20 )
    ggsave(paste0(program_name,"/boxplot_",i,"_wilcox_progeny.pdf"),width = 2,height = 3,family="serif")
  }
  
  ##
}

#UMAP----
##load data #all exp meta----
#exp<-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\1.Quality_Control\\raw_exp.csv",row.names = 1)
#meta<-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\1.Quality_Control\\overall_annotation.csv",row.names = 1)
#meta$id <-rownames(meta)#ID
#choose program data
exp<-exp[,rownames(meta)]

if(T){
  ##
  ##UMAP
  #读取为seurat 
  options(scipen = 100)#小数而非科学计数法
  library(Seurat)
  sce=CreateSeuratObject(counts = exp,project = "DSP",assay = "RNA")
  #sce@assays$RNA@counts[1:5,1:6]
  #add meta
  sce <- AddMetaData(sce,metadata = meta )
  library(RColorBrewer);library(ggpubr)
  colors <- colorRampPalette( brewer.pal(12, 'Paired') )( 20)
  
  
  ##plot-VlnPlot----
  
  for (i in c("TMA.position","Group","Type") ){
    print(i)
    print(
      VlnPlot(sce,group.by = i, #样本
              features = colnames(meta[13:16]), #纵坐标
              ncol = 4,log = T,pt.size =0,cols = colors)&  
        #theme_bw()&  
        theme(axis.title.x = element_blank(),        
              axis.text.x = element_text(color = 'black', size = 12),      #face = "bold",  
              #axis.text.y = element_text(color = 'black', face = "bold"),        
              axis.title.y = element_text(color = 'black',  size = 15),        
              panel.grid.major = element_blank(),        panel.grid.minor = element_blank(),        
              #panel.border = element_rect(color="black",size = 1.2, linetype="solid"),        
              panel.spacing = unit(0.12, "cm"),        
              plot.title = element_text(hjust = 0.5, face = "plain"),   #bold.italic     
              legend.position = 'none')
    )
    ggsave(paste0(program_name,"/VlnPlot_",i,".pdf"),width = 12,height = 3,family="serif")
  }
  
  ##SCP包-未安装
  #library(SCP)
  
  ##plot-UMAP----
  #DefaultAssay(sce) <- "RNA"
  sce_all <- sce %>%  NormalizeData() %>%
    FindVariableFeatures()  %>% 
    ScaleData() %>% RunPCA(npcs = 10) %>% #dims = 1:5
    #FindNeighbors(dims = 1:15) %>%
    #FindClusters( ) %>%#多则慢 #default 0.8 resolution = c(0.4, 0.6, 0.8, 1.0, 1.4)
    #RunTSNE(dims = 1:15) %>%
    RunUMAP(dims = 1:10)#慢
  #30个样本报错
  # Centering and scaling data matrix
  # |=================================================================================| 100%
  # Error in irlba(A = t(x = object), nv = npcs, ...) : 
  #   max(nu, nv) must be strictly less than min(nrow(A), ncol(A))
  #1:10 失败;PCA npcs = 10 成功
  
  for (i in c("TMA.position","Group","Type") ){
    print(i)
    print(
      DimPlot(sce_all, reduction = "umap",raster=T,group.by = i,pt.size = 5,
              cols = c("pink","skyblue","#46c1be",colors)
      )
      
      #labs(title = "before integration") 
    )
    ggsave(paste0(program_name,"/Umap_",i,".pdf"),width = 4,height = 3.5,family="serif")
  }
  
  ##
}

#转录因子-拟时序-hdWGCNA