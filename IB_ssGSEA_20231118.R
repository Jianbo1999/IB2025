
#20231118-IB-SSGSEA

rm(list=ls() )
setwd("G:\\IB\\231-202306\\transcriptomics\\ssGSEA")
#load("G:/IB/231-202306/transcriptomics/ssGSEA/IB_ssGSEA_20231118.RData")

#0. 加载数据----
exp <-read.csv("./data/All_gene_fpkm.csv")  #表达FPKM
exp <-exp[,-1]
library(limma)
exp_aver <- aggregate(.~Symbol, mean, data = exp) #重复基因取均值
rownames(exp_aver)<-exp_aver$Symbol;exp_aver<-exp_aver[,-1]

library(openxlsx)
ferroptosis <- read.xlsx("./data/markers.xlsx")  #铁死亡基因-转为list
ferroptosis_set<-list()
for (i in unique(ferroptosis$Pathway)){
  print(i);print(length(ferroptosis[ferroptosis$Pathway== i,]$Gene))
  ferroptosis_set[i]<-list(ferroptosis[ferroptosis$Pathway== i,]$Gene)
}
names(ferroptosis_set)[3]<-"UFA biosynthesis" #Unsaturated fatty acids biosynthetic process  UFA biosynthesis
names(ferroptosis_set)[4]<-"ROS" #Hallmark reactive oxygen species pathway ROS
#基因数太多？ 热图拼ssGSEA?/拼

mitogenes <- readxl::read_xls("./data/Human.MitoCarta3.0.xls",sheet = 4)   #线粒体基因
mitogenes_set <- split(mitogenes$Genes, mitogenes$MitoPathway)

#http://www.zhounan.org/ferrdb/current/operations/download.html 铁死亡驱动或者抑制基因
ferroptosis_suppressor<-read.csv("./data/ferroptosis_suppressor.csv") #238
ferroptosis_driver<-read.csv("./data/ferroptosis_driver_ferrdbv2.csv") #264


#1.0 hallmark 50 可视化----
#G:\bioinfo\github\FigureYa318GenesetDEDotplot(2) 
library(clusterProfiler)
geneset <- read.gmt("./data/h.all.v7.5.1.symbols.gmt")
geneset <- split(geneset$gene, geneset$term)
genesetInfo <- read.delim("./data/GenesetInfo.txt", sep = ",")
genesetInfo <- subset(genesetInfo, Classification != "")
levels = c("Immune", "Metabolism", "Signaling", "Proliferation")
genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
genesetInfo$Pathway <-tolower(gsub("HALLMARK_", "", levels(genesetInfo$geneset)))

##1.1 打分----
hallmark50_score <- gsva(as.matrix(exp_aver), #
                         geneset, #基因集
                    method = "ssgsea", kcdf = "Poisson", min.sz = 10)#


##1.2差异通路----
library(limma)
#library(edgeR)
#分组信息
group <-data.frame(name=c("C1","C2","C3","D1","D2","D3"),
                  group=rep(c('C', 'D'), each = 3, length.out = 6) )
#group_list <-c("C1","C2","C3","D1","D2","D3")#rep(c('C', 'D'), each = 3, length.out = 6) #A vs B 组 谁是给药组？  B vs A  #c(rep("A",3),rep("B",3))
#dgelist <- DGEList(counts = a549_exp, group = group)

exprSet <- hallmark50_score #log2
contra <-c('D-C')


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

##1.3 绘图----
# 准备画图数据
DPs <- read.table("./out/output_hallmark.txt", header = T)
plot.data <- DPs
plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
plot.data$Celltypes <- "IB vs Control"#factor(plot.data$Celltypes)
plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset) # 通路顺序与genesetInfo一致
plot.data$FDR<- cut(plot.data$adj.P.Val, breaks = c(0, 0.05, 0.5, 1),
                    include.lowest = T)
plot.data$FDR <- factor(as.character(plot.data$FDR),
                        levels = rev(levels(plot.data$FDR)))
levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))
##1.3.1 柱状图hallmark----
#类似NES:https://zhuanlan.zhihu.com/p/638107759

#-匹配factor
plot.data_df1 <-plot.data
plot.data_df1$Pathway<-factor(plot.data_df1$Pathway,levels = c(df1$y) )
plot.data_df1<-merge(plot.data_df1,genesetInfo,by="Pathway")
library(dittoSeq)
p1 <-ggplot(plot.data_df1,aes(-log10(adj.P.Val),Pathway))+geom_vline(xintercept = c(1.30103), linetype=2, linewidth=0.25)+  #1.30103
  geom_bar(stat = "identity", aes(fill=Classification ))+ #logFC
  geom_text(aes(label=round(logFC,2), x=-log10(adj.P.Val)+0.2 ),size=3)+  #logFC: log2 fold change limma 2^(logFC)
  labs(y='',x='-log10(adj.pvalue)')+ theme_classic()+#theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour = 'black', size = 10),
        plot.margin = margin(0,0,0,-0.05, "cm"))+
  scale_fill_manual(values = c("#0072B2","#F0E442","#009E73","#56B4E9","#E69F00"))
  #scale_fill_manual(values = dittoColors())
p1
p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) +
  geom_tile() +
  theme_classic() +
  #scale_fill_manual(values = class.color) +
  theme(axis.text = element_blank(), axis.title = element_blank(),legend.position = "right",
        axis.ticks = element_blank(), axis.line = element_blank())+
  coord_flip()+scale_fill_manual(values = c("#0072B2","#F0E442","#009E73","#56B4E9","#E69F00"))
p2
#p1+p2#p2 line 157
library(deeptime)
pdf( "./out/hallmark_bar.pdf", width = 7, height = 6,family="serif")
ggarrange2(p1, p2,nrow = 1,widths =c(2,0.05))
dev.off()

##1.3.2热图----
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

p3 #热图score
score_df<-as.data.frame(hallmark50_score)
score_df$Pathway <- factor(rownames(score_df), levels = genesetInfo$geneset) # 通路顺序与genesetInfo一致
levels(score_df$Pathway) <- tolower(gsub("HALLMARK_", "", levels(score_df$Pathway)))
score_df<-na.omit(score_df)  #删除NA通路  45

score_df[,1:6] <-scale(score_df[,1:6])
rownames(score_df) <-score_df$Pathway
 #排序

df_p <- t(score_df[,1:6])#t(scale( )) #Z-score
library (reshape2);library(dplyr)#数据转换
df_p <-df_p %>% melt() #长列表
colnames(df_p )<-c("sample","pathway","value")#sample pathway value

ggplot(df_p,aes(x=sample,y=pathway,fill=value,group=sample))+
  geom_tile(color="gray")+labs(x="",y="",fill="Score")+
  theme(legend.key.size = unit(0.15, "inches"),legend.title = element_text(size = 10),
        #axis.text.x =element_blank(),
        axis.ticks.x = element_blank(),panel.background = element_blank(),panel.grid = element_blank() )+ # + theme_bw()
  #scale_fill_gradientn(values = seq(0,1,0.2),colours = c('#6699CC','#FFFF99','#CC3333'))
  #scale_fill_viridis_c(option = "D",alpha=0.9)   #+facet_grid(~sample, scales = "free") #一模一样？组见差距较小
  scale_fill_distiller(palette = "Spectral") #scale_fill_gradientn(colors = brewer.pal(3, "Blues"))#scale_fill_brewer(palette="Set1")

library(RColorBrewer);display.brewer.all(type = "all")
display.brewer.all(type = "seq")

library(pheatmap)
pheatmap(score_df[,1:6])
hallmark_p <-pheatmap(score_df[,1:6],cluster_rows = F,cluster_cols  = F, show_colnames = T,
         border_color = "white",#cellwidth=10,cellheight = 10,
         gaps_col= 3, #table(exp_marker1$response)
         color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)) ,#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
         annotation_col=annotation_col,annotation_row = anno_row
)
hallmark_p
save_pheatmap_pdf(hallmark_p, "./out/hallmark_pheatmap.pdf",6,6) 

annotation_col #line 370
anno_row <- data.frame(
  Type=genesetInfo$Classification) #, size = 20, replace = TRUE  adj.P.Val
rownames(anno_row) <- tolower(gsub("HALLMARK_", "", levels(genesetInfo$geneset)))

#热图各组变化不明显！-放弃20231129

#1.4 il6-jak2-stat3  热图----
#1.4.1 提取基因查看热图
names(geneset[7])#87个
geneset[7][1]

#scale 样本
genes<-data.frame(unlist(geneset[7] ) )[,1]
exp_genes <-exp_aver[rownames(exp_aver) %in% genes,]
exp_genes <-exp_genes[which(rowSums(exp_genes) > 1),] #删除0行,和极小值行 57
exp_genes <-t(exp_genes)
exp_genes<-as.data.frame(scale(exp_genes) )

library("gplots");library(viridis)
hallmark_g <- pheatmap(t(exp_genes), #exp_genes[rownames(exp_genes) !=c("CD44","CD9"),],
          cluster_rows = T,cluster_cols  = F,gaps_col= 3,border_color = NA,
          #cluster_rows = F,cluster_cols  = T,gaps_row = 3,border_color = NA,
          #color= rev(colorRampPalette(brewer.pal(10, "Spectral"))(15) ),
          annotation_row = annotation_row_p,annotation_names_row = F,annotation_colors = ann_colors_p,annotation_col=annotation_col,annotation_names_col = F,
          #color=viridis(20, option = "D")[3:20], #G
          color=bluered(50)#library("gplots")
)
hallmark_g
save_pheatmap_pdf(hallmark_g, "./out/hallmark_IL6_JAK_STAT3_SIGNALING_scale_sample-anno-1.pdf",5,9) 

# 显著性注释
exp_noscale <-exp_aver[rownames(exp_aver) %in% genes,] #未scale
exp_noscale <-exp_noscale[which(rowSums(exp_noscale) > 1),] #删除0行,和极小值行

pvalue_g<-data.frame() #T.test
for (i in 1:nrow(exp_noscale) ){
  print(i)
  pwilcox <- t.test(as.numeric(exp_noscale[i,4:6] ),as.numeric(exp_noscale[i,1:3]) )
  fc <- mean(as.numeric( exp_noscale[i,4:6] ) )/mean( as.numeric(exp_noscale[i,1:3]) )#drug/control
  pvalue_g<-rbind(pvalue_g,data.frame(pvalue=pwilcox$p.value,FoldChange=fc,gene=rownames(exp_noscale)[i]) )
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
#"#FDBF6F" "#FF7F00""#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
     


#1.5 ##差异基因火山图中注释 ferroptosis genes----
if(T){
  library(ggplot2)
  library(ggsci)
  library(tidyverse)
  library(tidydr)
}



#load data
# 差异表格数据读入
Diff <- readxl::read_xls("G:/IB/231-202306/transcriptomics/Cmap/data/DEG_final-1.xls")[,c(2,15:17)]
colnames(Diff);colnames(Diff)[2:3]<-c("pvalue","log2FoldChange")
rownames(Diff)<-Diff$Symbol
#baseMean	log2FoldChange	lfcSE	pvalue	padj #公司数据未给padj-重新差异分析？



#1.5.1 ppi蛋白交互数据
ppiLink <- read.table("G:/IB/231-202306/transcriptomics/ssGSEA/29.火山图+ppi/9606.protein.links.v12.0.txt", header = T, sep = " ")
ppiAnn_fristP <- read.table("G:/IB/231-202306/transcriptomics/ssGSEA/29.火山图+ppi/9606.protein.info.v12.0.txt", header = F, sep = "\t")[, c(1,2)]
colnames(ppiAnn_fristP) <- c("protein1", "geneName1")
ppiAnn_secondP <- read.table("G:/IB/231-202306/transcriptomics/ssGSEA/29.火山图+ppi/9606.protein.info.v12.0.txt", header = F, sep = "\t")[, c(1,2)]
colnames(ppiAnn_secondP) <- c("protein2", "geneName2")
ppi_all <- merge(ppiLink, ppiAnn_fristP, by = "protein1") %>% merge(ppiAnn_secondP, by = "protein2") #巨大
ppi_all <- ppi_all[,c(4, 5, 3)]
head(ppi_all)



#火山图差异标记颜色标准设置
adjPval <- 0.05         # 矫正p值
aflogFC <- 1             # logFC
Significant <- ifelse((Diff$pvalue <adjPval & abs(Diff$log2FoldChange)>aflogFC), ifelse(Diff$log2FoldChange>aflogFC,"Up","Down"), "Not")
# 提取在Diff中出现的基因
ppi_line <- ppi_all[which(ppi_all$geneName1%in%rownames(Diff[!Significant=="Not",]) & ppi_all$geneName2%in%rownames(Diff[!Significant=="Not",])), ]


# 循环提取交互信息
line_res <- sapply(as.list(1:nrow(ppi_line)), function(x, Diff, ppi_line){
  index1 <- which(rownames(Diff)==ppi_line$geneName1[x])
  index2 <- which(rownames(Diff)==ppi_line$geneName2[x])
  return(c(Diff$log2FoldChange[index1], -log10(Diff$pvalue)[index1], Diff$log2FoldChange[index2], -log10(Diff$pvalue)[index2]))
}, Diff=Diff[!Significant=="Not",], ppi_line=ppi_line)
# 计算DC值：什么含义？
ppi_all <- ppi_all[ppi_all$geneName1 %in% rownames(Diff)[!Significant=="Not"],]
DC <- as.data.frame(table(ppi_all$geneName1))
rownames(DC) <- DC$Var1
DCnum <- ifelse(rownames(Diff) %in% rownames(DC), 1, 0)
DCnum[which(rownames(Diff) %in% rownames(DC))] <- (DC[rownames(Diff)[which(rownames(Diff) %in% rownames(DC))], ])[,2]
Diff$DCscore <- DCnum
Diff$DCscore <- (Diff$DCscore/max(Diff$DCscore))*3
#write.csv(Diff, "./out/Diff_result.csv", row.names = T)   # 将每个显著差异基因的dc值储存
# 根据DC值标记名称/ 标记自己感兴趣基因？
for_label <- Diff[!Significant=="Not",] #padj 
for_label <- for_label[abs(for_label$log2FoldChange) > 1.5 & for_label$pvalue < 0.05,] # 先根据logFC和padj进一步筛选 1e-64
for_label <- for_label[order(for_label$DCscore, decreasing = T)[1:20],]  # 剩余基因中前20个dc最大的基因
for_label$hubGene <- rownames(for_label)
Diff$hubGene <- rownames(Diff)
rm(ppiLink, ppiAnn_fristP, ppiAnn_secondP)

# 进行绘图
rownames(line_res) <- c("x", "y", "xend", "yend")
line_res <- t(line_res)%>%as.data.frame()
line_res$combined_score <- (ppi_line$combined_score/max(ppi_line$combined_score))*0.5
#开始绘制
p  <-  ggplot(Diff, aes(log2FoldChange, -log10(pvalue)))#pvalue  padj
p <- p+geom_segment(data=line_res, aes(x=x, y=y, xend=xend, yend=yend), color="green3", size=0.2, alpha=0.1)
p <- p+
  geom_point(aes(col=Significant, size=DCscore))+
  #scale_color_manual(values=c(ggsci::pal_nejm()(2)[2], "#838B8B", ggsci::pal_nejm()(2)[1]))+
  scale_color_manual(values=c("skyblue","grey","pink") )+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(adjPval)), colour="black", linetype="twodash",size=0.5)+
  geom_vline(aes(xintercept=aflogFC), colour="black", linetype="twodash",size=0.3)+
  geom_vline(aes(xintercept=-aflogFC), colour="black", linetype="twodash",size=0.3)+
  theme_dr(xlength = 0.1, ylength = 0.1) + theme(panel.grid=element_blank())+
  scale_size(range = c(0.5, 5))+
  ggrepel::geom_text_repel(
    aes(label = hubGene),
    data = for_label,
    color="black",max.overlaps=20,
    #label.size =0.01, 
    #max.overlaps = Inf,
     segment.size=0.4, #box.padding = 0.4,
    force = 10
  )
p #theme volid...
ggsave("./out/Diff_volcanoPlot-4.pdf", width = 6, height = 5)    # 保存文件时需要去除井号键

#1.6 抗原呈递相关基因变化？热图----
#基因集来源：

#2.0 mitogenes_set----

#3.0 ferroptosis_set----
library(GSVA)
ferro_score <- gsva(as.matrix(exp_aver), #
                         ferroptosis_set, #基因集
                         method = "ssgsea", kcdf = "Poisson", min.sz = 1)#
library(pheatmap);library(RColorBrewer)
pheatmap(ferro_score,#scale(ferro_score),
         color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(25) ),
         cluster_rows = F,cluster_cols  = F,gaps_col= 3)
##3.1 EDG----
library(limma)
fit <- lmFit(ferro_score,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$Pathway <-rownames(DEG)
# 输出差异分析表格
write.table(DEG, file = "./out/output_hallmark_ferro_score.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)
#差异结果不显著
##3.2 plot dot 
plot.data1 <-DEG
plot.data1$Celltypes <- "IB vs Control"#factor(plot.data$Celltypes)
plot.data1$FDR<- cut(plot.data1$adj.P.Val, breaks = c(0, 0.05, 0.5, 1),
                    include.lowest = T)
plot.data1$FDR <- factor(as.character(plot.data1$FDR),
                        levels = rev(levels(plot.data1$FDR)))

library(ggplot2)
ggplot(plot.data1, aes(x = Pathway, y = Celltypes, color = logFC, size = FDR#adj.P.Val#
)) +
  geom_point() +
  scale_color_gradient2(low = color[1], mid = color[2], high = color[3]) +
  theme_classic() +
  theme(axis.title.y = element_blank(),axis.title = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 14),axis.line = element_blank())+ #, legend.position = "top"
  coord_flip()
ggsave("./out/ferro_score_limma.pdf", width = 4, height = 5,family="serif")

###分开花热图和点图，合并顺序？
ferro_score1 <-as.data.frame(ferro_score)
ferro_score1$Pathway <-rownames(ferro_score1)
ferro_score1<-merge(ferro_score1,plot.data1,by="Pathway")
ferro_score1$FDR<- cut(ferro_score1$adj.P.Val, breaks = c(0, 0.05, 0.5, 1),
                     include.lowest = T)
ferro_score1$FDR <- factor(as.character(ferro_score1$FDR),
                         levels = rev(levels(ferro_score1$FDR)))

ferro_score1 <-ferro_score1[order(ferro_score1$FDR,decreasing = T),]
rownames(ferro_score1)<-ferro_score1$Pathway;ferro_score1<-ferro_score1[,-1]
ferro_score1$Pathway <-rownames(ferro_score1)
ferro_score1$Pathway <-factor(ferro_score1$Pathway,levels=rev(levels(ferro_score1$Pathway)) )

p3 <- ggplot(ferro_score1, aes(x = Pathway, y = Celltypes, color = logFC, size = FDR#adj.P.Val#
)) +geom_point() +
  #scale_colour_viridis_c(option = "D",direction = -1)+ 
  scale_color_gradient2(low = color[1], mid = color[2], high = color[3]) +
  theme_classic() +theme(axis.title.y = element_blank(),axis.title = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 14),axis.line = element_blank(), legend.position = "left")+ #, legend.position = "top"
  coord_flip()
p3
ggsave("./out/ferro_score_limma_p3-1.pdf", width = 4, height = 5,family="serif")

library(viridis)
p4 <-pheatmap(scale(ferro_score1[,1:6]),  #scale(ferro_score),
         #color= rev(colorRampPalette(brewer.pal(10, "Spectral"))(15) ),
         #color=viridis(20, option = "D")[3:20], #G
         color=bluered(25),#library("gplots")
         cluster_rows = F,cluster_cols  = F,gaps_col= 3,border_color = 'NA',
         annotation_row=annotation_row,annotation_col=annotation_col,
         annotation_colors=ann_colors
         ) 
p4
save_pheatmap_pdf(p4, "./out/ferro_score_limma_p4-3.pdf",6,5) #10,7


p5 <-pheatmap(scale(ferro_score1[,1:6]),  #scale(ferro_score),
              #color= rev(colorRampPalette(brewer.pal(10, "Spectral"))(10) ),
              #color=bluered(25),#library("gplots")
              color=viridis(20, option = "D")[5:20],
              cluster_rows = F,cluster_cols  = F,gaps_col= 3,border_color = 'NA')
p5
save_pheatmap_pdf(p5, "./out/ferro_score_limma_p5-1.pdf",4,4)

#注释信息

annotation_row <- data.frame(
  FDR=ferro_score1$FDR,logFC=ferro_score1$logFC) #, size = 20, replace = TRUE  adj.P.Val
rownames(annotation_row) <- rownames(ferro_score1)
annotation_col <-data.frame(
  Group=c(rep("Control",3),rep("Drug",3)) )
rownames(annotation_col) <-colnames(ferro_score1[,1:6])
ann_colors=list(FDR=c('[0,0.05]'='#FDDCA9','(0.05,0.5]'="#DDC9E3"),
                Group=c("Control"="skyblue","Drug"="pink" ))

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

#alpha("purple",0.5) "#A020F0B2"  '#A020F080''
#'#1B9E77'  '#D95F02'
#"#7294D4" "#FB9A99FF" 
#'#AA2E2F'   '#7093E8'    '#F4A6AB' '#DDC9E3'  '#F98E1F'   '#A03230'

##3.2 pheatmap: all ferroptosis geneset----
all_ferro_exp <-exp_aver[rownames(exp_aver) %in% ferroptosis$Gene ,]#268-6

#all_ferro_exp<-as.data.frame(t(all_ferro_exp))#scale基因在样本中
#all_ferro_exp <-as.data.frame(scale(all_ferro_exp) );boxplot(all_ferro_exp)
  #2~-2
# library(dplyr)
# all_ferro_exp <-all_ferro_exp |> 
#   mutate(across(where(is.numeric),
#                 ~scales::rescale(.x,
#                                  to=c(-2,2)))) #归一化-2~2
#计算差异wilcox

pvalue<-data.frame()
for (i in 1:nrow(all_ferro_exp) ){
  print(i)
  pwilcox <- wilcox.test(as.numeric(all_ferro_exp[i,1:3] ),as.numeric(all_ferro_exp[i,4:6]) )
  pvalue<-rbind(pvalue,pwilcox$p.value)
} #全部无显著性-wilcox
#all_ferro_exp$pvalue <-pvalue[,1]

#提取某基因boxplot ----
gene <- "CD274"  #ARPPARA ATF3 NFE2L2/nrf2 TP53 ATF4 AIFM2/FSP1 FH CASP1  NLRP3 GSDME"GAPDH" FTH1 FTL STAT3 NCOA4 GPX4 SLC7A11 SLC47A1 PPARA
#CD274增加，不符合趋势
if(T){
  df_gene <-exp_aver[rownames(exp_aver) %in% gene,]
  df_gene <-as.data.frame(t(df_gene))
  df_gene <-as.data.frame(scale(df_gene) )  #scale样本
  df_gene$Group <- c(rep("Control",3),rep("Drug",3))
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
    stat_compare_means(method = "t.test",size=5,show.legend= F,label = "p.format",label.x =0.75,label.y =max(df_gene$Gene)-0.06)#+ #,label.y =max(df_gene$Gene)-0.55,label.x =1.5,label.y = 1.3  0.8 ,label.y = 4 #非参数检验kruskal.test not anova
    ylim(NA,max(df_gene$Gene)+0.5)
    ggsave(paste0("./out/boxplot_Ttest_",gene,".pdf"),width = 2,height = 3,family="serif")
  
}
ggsave(paste0("./20240131out/boxplot_",gene,"_Ttest.pdf"),width = 2,height = 3,family="serif")



p6 <-pheatmap(scale(all_ferro_exp),  #scale(ferro_score),
         #color= rev(colorRampPalette(brewer.pal(10, "Spectral"))(10) ),
         #color=bluered(25),#library("gplots")
         #color=viridis(20, option = "D")[5:20],
         fontsize = 10, 
         fontsize_row = 8,  
         fontsize_col = 10,#scale = "row",#column
         cluster_rows = T,cluster_cols  = F,gaps_col= 3,border_color = 'NA')
p6
save_pheatmap_pdf(p6, "./out/all_ferro_exp-2.pdf",4,24)



#! Iron take
table(ferroptosis$Pathway)#Iron intake 17


##3.3 GSEAvis----
#自定义GSEA https://mp.weixin.qq.com/s/wFY2zo01FYfN0h0K2REXug

#4.0 重新差异分析20240131----
#4.0 limma----
## 导入R包
library(limma)
library(dplyr)
#https://zhuanlan.zhihu.com/p/437622149?utm_id=0
#https://blog.csdn.net/weixin_45161743/article/details/103536340
#样本注释
list <- c(rep("Control", 3), rep("Drug",3)) %>% factor(., levels = c("Control", "Drug"), ordered = F)
head(list)
list <- model.matrix(~factor(list)+0)  #把group设置成一个model matrix
colnames(list) <- c("Drug", "Control")
df.fit <- lmFit(exp_aver, list)  ## 数据与list进行匹配

dge <- DGEList(counts = exp_aver)#library(edgeR)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE) #会自动计算log(cpm)值
#拟合线性模型
fit <- lmFit(v, design)
#针对给定的对比计算估计系数和标准误差
fit2 <- contrasts.fit(fit, df.matrix)
#计算出t统计量，F统计量和差异表达倍数的对数
fit2 <- eBayes(fit2)
#提取差异基因
allDEG <- topTable(fit2,n = Inf, adjust = "fdr")
allDEG <- na.omit(allDEG)
# padj = 0.05
# foldChange= 1
# diff_signif = allDEG[(allDEG$adj.P.Val < padj & abs(allDEG$logFC)>foldChange),]                    
# diff_signif = diff_signif[order(diff_signif$logFC),] #adj.P.Val   180 +newgene
# save(diff_signif, file = './20210131out/limma_diff.Rdata')

allDEG$Significant <- ifelse((allDEG$adj.P.Val <adjPval & abs(allDEG$logFC)>aflogFC ), ifelse(allDEG$logFC>aflogFC,"Up","Down"), "Not")
allDEG$Symbol <-rownames(allDEG) #删除新基因
allDEG$Symbol <-gsub("_.*", "", allDEG$Symbol   )#1505
allDEG <- subset(allDEG ,allDEG$Symbol != "Homo")#19957

diff_signif<-subset(allDEG,allDEG$Significant !="Not")#151
table(diff_signif$Significant)#Down   Up  37  114


#4.1 VOLCANO ----
# Diff_1 <-Diff  
# Diff_1$Symbol <- gsub("_.*", "", Diff_1$Symbol  )#1505
# Diff_1 <- subset(Diff_1 ,Diff_1$Symbol != "Homo")#1426
# Diff_1$pvalue  <- ifelse(Diff_1$pvalue <= 1.430066e-50,1.430066e-50,Diff_1$pvalue) #pvalue
# #提取铁死亡相关差异基因

Diff_1 <-allDEG;colnames(Diff_1)
colnames(Diff_1)[1]<-"log2FoldChange"
colnames(Diff_1)[4]<-"pvalue"

#标记
#for_label_2 <- subset(Diff_1[Diff_1$adj.P.Val <0.05,],Diff_1[Diff_1$adj.P.Val <0.05,]$Symbol %in% ferroptosis_driver$symbol)#ferroptosis$Gene
#for_label_2 <-  subset(for_label_2,for_label_2$Significant !="Not")#显著5
for_label_2 <- subset(Diff_1,Diff_1$Ferro != "None")#ferroptosis$Gene

#shape标记铁死亡驱动/抑制基因
#Diff_1<-Diff_1[,-9]
#Diff_1$Ferro <-ifelse(Diff_1$Symbol %in% ferroptosis_driver$symbol,"Driver", "None")
#Diff_1$Ferro <-ifelse(Diff_1$Symbol %in% ferroptosis_suppressor$symbol,"Suppressor", "None")
Diff_1$Ferro <-ifelse(Diff_1$Symbol %in% ferroptosis_driver$symbol,"Driver", 
                      ifelse(Diff_1$Symbol %in% ferroptosis_suppressor$symbol,"Suppressor", "None")
                      )

library(ggplot2)
#pdf("./20240131out/Diff_volcanoPlot-_p-2.pdf", width = 7, height = 7,family="serif")
ggplot(Diff_1, aes(log2FoldChange, -log10(adj.P.Val), shape=Ferro  ) )+
  geom_point(aes(col=Significant  ),size=3,alpha=0.8)+ #, size= abs(log2FoldChange)
  #scale_color_manual(values=c(ggsci::pal_nejm()(2)[2], "#838B8B", ggsci::pal_nejm()(2)[1]))+
  scale_color_manual(values=c("skyblue","grey","pink") )+
  labs(title = " ",x="log2(FoldChange)",y="-log10(adj.pvalue)" )+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),legend.position = "top")+
  geom_hline(aes(yintercept=-log10(adjPval)), colour="black", linetype=2,size=0.5)+
  geom_vline(aes(xintercept=aflogFC), colour="black", linetype=2,size=0.3)+
  geom_vline(aes(xintercept=-aflogFC), colour="black", linetype=2,size=0.3)+
  theme_classic(base_size = 14)+#ylim(0,4.5)+

  ggrepel::geom_text_repel(
    aes(label = Symbol),
    data = for_label_2,
    color="black",max.overlaps=80,
    #label.size =0.01, 
    #max.overlaps = Inf,
    segment.size=0.3, box.padding = 0.6#,force = 10
  )
#dev.off()
#中间有缝-重新差异分析！
table(Diff_1$Significant)#Down  Not   Up  37 19806  114
#ggsave("./out/Diff_volcanoPlot-ferro-1.pdf", width = 8, height = 6,family="serif")    # 保存文件时需要去除井号键
ggsave("./20240131out/Diff_volcanoPlot-_p-1.pdf", width = 7, height = 7,family="serif")    # 保存文件时需要去除井号键
#suppressor  driver slc7a11均有

#20240222-up、down各显著driver/supposer
library(ferroviz);data(ferrdbhs)
Diff_2 <-Diff_1
Diff_2$FoldChange <-2^(Diff_2$log2FoldChange)
plot(density(Diff_2$FoldChange))
plot(density(Diff_2$log2FoldChange))
Diff_2$label <-ifelse(Diff_2$adj.P.Val <0.05 & Diff_2$Symbol %in% ferrdbhs[ferrdbhs$class=="driver",]$hs.gene,Diff_2$Symbol,
                     ifelse(Diff_2$adj.P.Val <0.05 & Diff_2$Symbol %in% ferrdbhs[ferrdbhs$class=="suppressor",]$hs.gene,Diff_2$Symbol,NA )
                     )#up down 太少
#Diff_2$Type <-ifelse(Diff_2$adj.P.Val <0.05 & Diff_2$Symbol %in% ferrdbhs[ferrdbhs$class=="suppressor",]$hs.gene,Diff_2$Symbol,NA )



for_label_3 <-subset(Diff_2,Diff_2$label != "NA") #76 #Diff_2[Diff_2$label != "NA",]
#up - driver; down -suppressor #8; FC>1 45
for_label_4 <-subset(for_label_3,
                     (for_label_3$log2FoldChange>0.5 & for_label_3$Ferro=="Driver") | (for_label_3$log2FoldChange< -0.5 & for_label_3$Ferro=="Suppressor") ) #76 #Diff_2[Diff_2$label != "NA",]

colnames(Diff_2)
ggplot(Diff_2, aes(log2FoldChange, -log10(adj.P.Val), shape=Ferro,color=Change  ) )+ #Significant
  geom_point(size=3,alpha=0.8)+ #aes(col=Significant  ),  , size= abs(log2FoldChange)
  #scale_color_manual(values=c("skyblue","grey","pink") )+
  scale_color_brewer(palette="Paired",direction = -1)+# display.brewer.pal(n=8,name="Paired")
  #ggsci::scale_color_npg()+
  labs(title = " ",x="log2(FoldChange)",y="-log10(adj.pvalue)" )+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),legend.position = "top")+
  geom_hline(aes(yintercept=-log10(adjPval)), colour="black", linetype=2,size=0.5)+
  geom_vline(aes(xintercept=0.5), colour="black", linetype=2,size=0.3)+ #aflogFC
  geom_vline(aes(xintercept=-0.5), colour="black", linetype=2,size=0.3)+
  theme_classic(base_size = 14)+
  ggrepel::geom_text_repel(
    aes(label = Symbol,color=Ferro), #,color=Ferro
    data = for_label_4,show.legend = F,
    #color="black",
    max.overlaps=80,
    #label.size =0.01, 
    #max.overlaps = Inf,
    segment.size=0.3, box.padding = 0.6 )#,force = 10
ggsave("./20240131out/Diff_volcanoPlot-ferro-20240222-ds-lo2fc0.5-up523-down683.pdf", width = 8, height = 6,family="serif")    # 保存文件时需要去除井号键

##重新阈值FC log2FoldChange 0.5
Diff_2$Change <- ifelse((Diff_2$adj.P.Val <adjPval & abs(Diff_2$log2FoldChange)>0.5), ifelse(Diff_2$log2FoldChange>0.5,"Up","Down"), "Not")
table(Diff_2$Change)#lo2FC down 523 up 683
table(Diff_2$Significant)#lo2FC 1  down 37 up 144


#4.2 KEGG enrichment----
#所有差异基因富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
                         #[diff_signif$Significant=="Down",]
input <- bitr(diff_signif$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")#1271

if(T){
  kegg_all <- enrichKEGG(gene = input$ENTREZID, #input$ENTREZID
                         organism ='hsa', #"hsa" #Homo sapiens (human)
                         pAdjustMethod = "BH",keyType = "kegg",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         use_internal_data =FALSE)
  kegg_all <- setReadable(kegg_all,keyType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
  kegg_all_res <- as.data.frame(kegg_all@result)
  kegg_all_res_p <- subset(kegg_all_res,kegg_all_res$pvalue <0.05)#57  铁死亡不显著！ down-7  up-55 pvaule
  
}
#write.csv(kegg_all_res,"./20240131out/kegg_all_res.csv")

##可视化APEAR----
#https://mp.weixin.qq.com/s/upJj5SPLaZHrtvl1qOSjTg
#install.packages("aPEAR") 
library(data.table);library(tidyverse) 
library(clusterProfiler) ;library(DOSE) ;library(org.Hs.eg.db);library(aPEAR)
data(geneList)


# 设置另一个随机数种子，用于后续的可视化过程
set.seed(654824)
# 创建富集分析的网络图，这里使用enrich@result作为输入数据
enrichmentNetwork(kegg_all_res_p)
enrichmentNetwork(kegg_all)

# 设置随机数种子，用于可视化的稳定性
set.seed(348934)
# 创建基于p值的富集网络图，这里指定了颜色类型为p值，pCutoff为p值的截断阈值
enrichmentNetwork(kegg_all@result, colorBy = 'pvalue', colorType = 'pval', pCutoff = -5)


##pathfindR 富集通路聚类、打分、比较与网络可视化----
#https://mp.weixin.qq.com/s/tpVAgy7wVwP-mjh6FFG_Ug
library(pathfindR) #安装失败
input_df<-allDEG[,c(8,1,5)]
output_df <- run_pathfindR(input_df)

#4.3 GSEA----

eg <- bitr(allDEG$Symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
allDEG1 <- allDEG[allDEG$Symbol %in% eg$SYMBOL,]
colnames(eg)[1] <-"Symbol"

allDEG1 <-merge(allDEG1,eg,by="Symbol")
genelist <-  allDEG1[allDEG1$ENTREZID,]$logFC#所有基因-18866
names(genelist) = allDEG1$ENTREZID
genelist = sort(genelist, decreasing = TRUE)

kkgsea <- gseKEGG(geneList     = genelist,
                  organism     = 'hsa', 
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 1,
                  pAdjustMethod = "none" ) #进行gseKEGG富集分析 
kkgsea=setReadable(kkgsea,keyType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
kkgsea_res <- kkgsea@result #331
kkgsea_res_p <- kkgsea_res[kkgsea_res$p.adjust <0.05,]#39  7-NES < 0
#write.csv(kkgsea_res,"./20240131out/gsea_res.csv")


library(enrichplot)
gseaplot2(kkgsea,geneSetID ="hsa04392",color = "red" ) #
gseaplot2(
  kkgsea, #gseaResult object，即GSEA结果
  geneSetID = "hsa04392",#1,#"hsa05168",#富集的ID编号
  title = "", #标题
  color = "pink",#GSEA线条颜色
  base_size = 11,#基础字体大小
  rel_heights = c(1.5, 0.5, 1),#副图的相对高度
  subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
  pvalue_table = T, #是否添加 pvalue table
  ES_geom = "line" )#running enrichment score用先还是用点ES_geom = "dot"

##4.3.2 铁死亡基因driver suppression 分开GSEA----
#ferroptosis_driver$symbol 基因集
geneset_ferro<-list(driver=ferroptosis_driver$symbol,suppressor=ferroptosis_suppressor$symbol)

genelist1 <-  allDEG1$logFC#所有基因-18866
names(genelist1) = allDEG1$Symbol
genelist1 = sort(genelist1, decreasing = TRUE)
genelist1 <- genelist1[genelist1 != 0]

geneset_ferro <- data.frame(term = c(rep("Ferroptosis driver",nrow(ferroptosis_driver) ),
                               rep("Ferroptosis suppressor",nrow(ferroptosis_suppressor) ) ),
                               
                      gene = c(ferroptosis_driver$symbol,ferroptosis_suppressor$symbol))
head(geneset_ferro)

library(clusterProfiler)
set.seed(123456)
my_GSEA <- GSEA(genelist1, TERM2GENE=geneset_ferro, #symbol 
                verbose=FALSE,nPerm = 10000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)

#plot
library(enrichplot);library(ggplot2)
my_GSEA@result[["ID"]]
gseaplot2(my_GSEA,geneSetID =c("Ferroptosis driver"),color = "red" ) #
gseaplot2(
  my_GSEA, #gseaResult object，即GSEA结果
  geneSetID = "Ferroptosis driver",#1,#"hsa05168",#富集的ID编号
  title = "", #标题
  color = "pink",#GSEA线条颜色
  base_size = 14,#基础字体大小
  rel_heights = c(1.5, 0.5, 1),#副图的相对高度
  subplots = 1:2, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图，subplots=1#只要第一个图
  pvalue_table = T, #是否添加 pvalue table
  ES_geom = "line" )#
ggsave("./20241224out/gseaplot_Ferroptosis_driver.pdf",family="serif",width = 5.5,height = 4)

#GSEAvis
library(GseaVis)
gseaNb(object = my_GSEA,
       geneSetID = "Ferroptosis driver", termWidth = 35,#hsa00310 hsa00480
       #newGsea = F,addPoint = F,
       #rmSegment = T, #移除红色标记线
       #rmHt = T,#移除热图
       addPval = T,#Add NES and Pvalue
       pvalX = 0.45,pvalY = 0.8, #调整标签位置和颜色
       pCol = 'red',pHjust = 0,subPlot = 2)
ggsave(paste0("./20241224out/GSEAvis_","Ferroptosis_driver",".pdf"),width = 4,height = 4)
ggsave(paste0("./20241224out/GSEAvis_","Ferroptosis_driver_1",".pdf"),width = 3.6,height = 3.7)
#new
gseaNb(object = my_GSEA,
       geneSetID = "Ferroptosis driver", termWidth = 35,#hsa00310 hsa00480
       #newGsea = F,addPoint = F,
       #rmSegment = T, #移除红色标记线
       #rmHt = T,#移除热图
       newGsea = T,#htHeight=0,
       addPval = T,#Add NES and Pvalue
       pvalX = 0.45,pvalY = 0.8, #调整标签位置和颜色
       pCol = 'red',pHjust = 0,subPlot = 2)

#添加表达热图
gseaNb(object = my_GSEA,
       geneSetID = "Ferroptosis driver", termWidth = 35,#hsa00310 hsa00480
       #newGsea = F,addPoint = F,
       #rmSegment = T, #移除红色标记线
       #rmHt = T,#移除热图
       #newGsea = T,
       add.geneExpHt = T,exp=exp,ght.geneText.size = 8,
       addPval = T,#Add NES and Pvalue
       pvalX = 0.45,pvalY = 0.8, #调整标签位置和颜色
       pCol = 'red',pHjust = 0,subPlot = 1)
ggsave(paste0("./20241224out/GSEAvis_","Ferroptosis_driver_pheatmap",".pdf"),width = 8,height = 5)

#5 兴趣基因热图----
gene_interest <-c("B2M","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB",
                  "HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1",
                  "HLA-DRA","HLA-DRB1","HLA-E","HLA-F","HLA-G","TAP1","TAP2","TAPBP" ) #"JAK2","STAT3","STAT4","CD274",

gene_interest <-c("JAK2","STAT3","STAT4","CD274","IL12RB1","IL12RB2")
#TCA
gene_interest <-c("PDHA1","PDHB","CS","ACO2","IDH3A","OGDH","DLST","DLD","SDHA","SDHB","SDHC","SDHD","SDHAF1","FH","MDH2")
#glycolysis
gene_interest <-c("HK2","HK1","GPI","PFK1","TPI","GAPDH","PGK","PGAM1","ENO","PKM1","PKM2","LDHA","LDHB",
                  "ENO1","ENO2","PFKFB3","PFKFB4")


if(T){
  exp_interest <-exp_aver[rownames(exp_aver) %in% gene_interest,]
  exp_interest <-exp_interest[which(rowSums(exp_interest) > 1),] #删除0行,和极小值行 57
  exp_interest <-t(exp_interest)
  exp_interest<-as.data.frame(scale(exp_interest) )
  
  library("gplots");library(viridis);library(pheatmap);library(RColorBrewer)
  
  exp_int <- exp_aver[rownames(exp_aver) %in% gene_interest,]
  exp_int <- exp_int[which(rowSums(exp_int) > 1),]
  
  pvalue_interest <-data.frame() #T.test
  # for (i in 1:nrow(exp_int) ){
  #   print(i)
  #   pwilcox <- t.test(as.numeric(exp_int[i,4:6] ),as.numeric(exp_int[i,1:3]) )
  #   fc <- mean(as.numeric( exp_int[i,4:6] ) )/mean( as.numeric(exp_int[i,1:3]) )#drug/control
  #   pvalue_interest<-rbind(pvalue_interest,data.frame(pvalue=pwilcox$p.value,FoldChange=fc,Symbol=rownames(exp_int)[i] ) )
  # } #wilcox.test 均不显著；t.test 较多显著性-3个, wilcox-0
  #对样本scale后计算,0则NA;有0值
  
  
  pvalue_interest<-Diff_1[Diff_1$Symbol %in% gene_interest,]
  pvalue_interest$FoldChange <- 2^pvalue_interest$log2FoldChange
  pvalue_interest$FC <-ifelse(2^pvalue_interest$log2FoldChange <1,"Decrease","Increase")
  
  
  #手动加上显著性？
  pvalue_interest$Sign. <- ifelse(pvalue_interest$pvalue>= 0.05,"ns",
                                  ifelse(pvalue_interest$pvalue < 0.01,"**","*") )
  # pvalue_interest$FC <-ifelse(pvalue_interest$FoldChange <=0.5,"Decrease",
  #                             ifelse(pvalue_interest$FoldChange >1,"Increase","Unchange") )
  rownames(pvalue_interest) <-pvalue_interest$Symbol#pvalue_interest$gene
  #合并差异分析结果！
  #pvalue_interest<-merge(pvalue_interest,Diff_1[Diff_1$Symbol %in% gene_interest,],by="Symbol")
  
  anno_row_p <- data.frame(
    Sign.=pvalue_interest$Sign.,FoldChange=pvalue_interest$FC) #, size = 20, replace = TRUE  adj.P.Val
  rownames(anno_row_p) <- rownames(pvalue_interest)
  ann_colors_p=list(FoldChange=c('Decrease'="skyblue",'Increase'="#4a6fe3"),  #'Unchange'='#CAB2D6',
                    Group=c("Control"="#7EC7A7","Drug"="#EDA065"),
                    Sign.=c("ns"="grey","*"="pink",'**'="#D33F6A" ))
  
  hallmark_interest <- pheatmap(t(exp_interest), #exp_genes[rownames(exp_genes) !=c("CD44","CD9"),],
                                cluster_rows = T,cluster_cols  = F,gaps_col= 3,border_color = NA ,
                                # #cluster_rows = F,cluster_cols  = T,gaps_row = 3,border_color = NA,
                                color= rev(colorRampPalette(brewer.pal(10, "Spectral"))(15) ),#PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 Spectral PRGn RdYlBu
                                #colorRampPalette(c("navy", "white", "firebrick3"))(50),#
                                annotation_names_row = F,annotation_row = anno_row_p,anno_names_row = F,annotation_colors = ann_colors_p,annotation_col=annotation_col,
                                # #color=viridis(20, option = "D")[3:20], #G
                                # color=bluered(50)#library("gplots")
  )
  hallmark_interest  #brewer.pal.info 
}
write.csv(t(exp_interest),"./out/TCA.csv")

#注释FC, pvalue
#save_pheatmap_pdf(hallmark_interest, "./20240131out/heatmap-HLAs-anno-1-2.pdf",5,4.5) 
#save_pheatmap_pdf(hallmark_interest, "./out/heatmap-srps-anno-1.pdf",5,2) #5,4.6
#save_pheatmap_pdf(hallmark_interest, "./20240131out/heatmap-tca-2.pdf",5,4.6) #5,4.6
#save_pheatmap_pdf(hallmark_interest, "./20240131out/heatmap-gly-gene.pdf",5,3) #5,4.6

#6.0 兴趣代谢物boxplot----
#代谢物boxplot
##6.0 加载合并离子模式非靶向表达谱
untarget_exp <-read.csv("G:\\IB\\231-202306\\Untargeted_metabolomics\\202312\\周建波_20231218\\处理后数据.csv")
untarget_exp$name <- gsub("w.o.MS2..","",untarget_exp$name )
untarget_exp<- aggregate( . ~ name,data=untarget_exp, mean)#合并重复值-平均
rownames(untarget_exp)<-untarget_exp$name;untarget_exp<-untarget_exp[,-1]
#untarget_exp<-untarget_exp[,-c(X0C,X2E)]#X0C  X2E 异常样本


meta_interest<-c("Glutathione..reduced.")#Glutathione..reduced. GSH #Glutathione..oxidized. GSSG # Alpha.Tocopherol CYSTEINE #Glutathione reduced ? oxidized?

if(T){
  df_gene <-untarget_exp[rownames(untarget_exp) %in% meta_interest,]
  df_gene <-as.data.frame(t(df_gene))
  df_gene <-as.data.frame(scale(df_gene) )  #scale样本
  df_gene$Group <- c(rep("Control",6),rep("Drug",6))
  colnames(df_gene)[1] <-"Gene"
  library(ggplot2);library(ggpubr)
  ggplot(df_gene,#
         aes(Group,Gene,color=Group))+  #age_group
    geom_boxplot(alpha=1,width=0.45,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    geom_point(alpha=0.7,size=1.5,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.8))+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    scale_color_manual(values =  c("#33CCFF","#FF66CC"))+ #"#EDA065","#66CCFF","#7EC7A7" blue c("#0066CC","#66CCFF","#9DC3E6")
    theme_bw() +labs(color="",title = meta_interest,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
    theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
          panel.grid = element_blank(),legend.key = element_blank(),legend.position="none" ) +
    stat_compare_means(method = "t.test",size=5,show.legend= F,label = "p.format",label.x =1,label.y =max(df_gene$Gene)-0.07
                       )#+ ylim(-1,1)
  
  ggsave(paste0("./20240131out/boxplot_",meta_interest,"_Ttest-0.pdf"),width = 2,height = 3,family="serif")
  
}
#c("#0066CC","#66CCFF","#9DC3E6")"#33CCFF" "#CCCCFF""#CC99FF" #00CC66
##6.1 兴趣代谢物
# meta_interest<-c("Citrate","Pyruvate","Acetyl coA","ketoglutarate","Succinate","Fumarate","Malate","Oxaloacetate")
# meta_interest<-c(meta_interest,toupper(meta_interest),tolower(meta_interest) )
# meta_exp <-untarget_exp[untarget_exp$name %in% meta_interest,]
# grepl(meta_interest[1],untarget_exp$name,ignore.case=TRUE)
# table(grepl(meta_interest[1],untarget_exp$name,ignore.case=TRUE) )
# meta_exp <-untarget_exp[grep(meta_interest[18],untarget_exp$name,ignore.case=TRUE,value=TRUE)==T,]
# grep(meta_interest[1],untarget_exp$name,ignore.case=TRUE,value=TRUE)
# 
# grepl(meta_interest[1],untarget_exp$name)
# table(grepl(meta_interest[1],untarget_exp$name) )

##直接在EXCEL中搜索-加载TCA数据----
untarget_exp_TCA <-read.csv("G:\\IB\\231-202306\\Untargeted_metabolomics\\202312\\周建波_20231218\\TCA.csv")
rownames(untarget_exp_TCA)<-untarget_exp_TCA$name
untarget_exp_TCA<-untarget_exp_TCA[,-1]
colnames(untarget_exp_TCA)<-c(paste0("C",1:6), paste0("D",1:6) )
write.csv(untarget_exp_TCA,"./out/TCA_matabolites.csv")
##绘图----
exp_int<- untarget_exp_TCA#t( scale(t(untarget_exp_TCA))  ) #z-score 样本标准化

#run
pvalue_interest <-data.frame() #T.test
for (i in 1:nrow(exp_int) ){
  print(i)
  pwilcox <- t.test(as.numeric(exp_int[i,7:12] ),as.numeric(exp_int[i,1:6]) )
  fc <- mean(as.numeric( exp_int[i,7:12] ) )/mean( as.numeric(exp_int[i,1:6]) )#drug/control
  pvalue_interest<-rbind(pvalue_interest,data.frame(pvalue=pwilcox$p.value,FoldChange=fc,Symbol=rownames(exp_int)[i] ) )
}

pvalue_interest$Sign. <- ifelse(pvalue_interest$pvalue>= 0.05,"ns",
                                ifelse(pvalue_interest$pvalue < 0.01,"**","*") )
pvalue_interest$FC <-ifelse(pvalue_interest$FoldChange <=2/3,"Decrease",
                            ifelse(pvalue_interest$FoldChange >=1,"Increase","Unchange") )
rownames(pvalue_interest) <-pvalue_interest$Symbol

anno_row_p <- data.frame(
  Sign.=pvalue_interest$Sign.,FoldChange=pvalue_interest$FC) #, size = 20, replace = TRUE  adj.P.Val
rownames(anno_row_p) <- rownames(pvalue_interest)
ann_colors_p=list(FoldChange=c('Unchange'='#CCECFF','Decrease'="skyblue",'Increase'="#4a6fe3"),  #'Unchange'='#CAB2D6',
                  Sign.=c("ns"="grey","*"="pink",'**'="#D33F6A" ) )

hallmark_interest <- pheatmap( t( scale(t(untarget_exp_TCA))  ), #t(exp_int) #此处scale,而前面不
                              cluster_rows = T,cluster_cols  = F,gaps_col= 6,border_color = NA ,
                              # #cluster_rows = F,cluster_cols  = T,gaps_row = 3,border_color = NA,
                              color= rev(colorRampPalette(brewer.pal(10, "Spectral"))(15) ), #PiYG 紫色绿色  RdYlGn Spectral PRGn RdYlBu
                              annotation_names_row = F,annotation_row = anno_row_p,anno_names_row = F,annotation_colors = ann_colors_p#,annotation_col=annotation_col,
                              # #color=viridis(20, option = "D")[3:20], #G
                              # color=bluered(50)#library("gplots")
)
hallmark_interest
##绘图结束

#绘图标准化，计算p FC 使用峰强度----save
#save_pheatmap_pdf(hallmark_interest, "./20240131out/heatmap-tca-metabo-2.pdf",5,2) #5,4.6


#7.0 ferroviz可视化----
##7.1安装----
library(devtools)
#install_github("cdesterke/ferroviz")
install_local("G:\\IB\\231-202306\\transcriptomics\\ssGSEA\\data\\ferroviz-main\\ferroviz-main")

#7.2 run----
#volcano
library(ferroviz)
data(resulths)
data(ferrdbhs)
resulths<- allDEG#limma结果"logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"  
hsvol(resulths,ferrdbhs,color="grey",size=10,x=0.15,label=1.2)+theme_bw()
#自己画ferrdbhs


#barplot
p2<-barploths(resulths,ferrdbhs,fc=0.25,size=16)
p2
#p2+theme_bw()
p2+guides(#color = FALSE, #删除color的legend
          fill = guide_legend(title="Type",#修改图例标题
                              title.theme = element_text( size = 15),
                              override.aes = list(colour = NA,alpha = 1,label = "")) #删除有颜色的legend ,
)+theme_bw()+theme(axis.text = element_text(size = 15),axis.title.y = element_text(size = 15),legend.text = element_text(size = 13),
                   legend.position = "top")
#ggsave("./20240131out/barplot-2-1.pdf",family="serif",width = 4,height = 5)#3.2 4

#+卡方检验
if(T){
  df_ferr <-p2$data  #resulths[(abs(resulths$logFC) >1 & resulths$P.Value <0.05) ,]#171  #p2$data# 
  table(df_ferr$ferroptosis)#38 46 #如何匹配
  
  # result <- df_ferr %>% 
  #   group_by(Significant, regulation) %>%
  #   summarize(count = n()) %>%
  #   mutate(proportion = count / sum(count) )
  # print(result)# 输出结果
  #https://blog.csdn.net/hx2024/article/details/134548265
  table(df_ferr$ferroptosis)
  table(df_ferr[df_ferr$regulation=="up",]$ferroptosis)[[1]]#27 26
  table(df_ferr[df_ferr$regulation=="down",]$ferroptosis)#11 20
  result <-data.frame(up=c(as.numeric(table(df_ferr[df_ferr$regulation=="up",]$ferroptosis))),down=c(as.numeric(table(df_ferr[df_ferr$regulation=="down",]$ferroptosis))))
  rownames(result)<-c("driver","suppressor")
  
  #执行检验
  #https://zhuanlan.zhihu.com/p/612852391?utm_id=0
  # 创建一个 2x2 的列联表 
  data <- matrix(c(27, 26,  11,20), nrow = 2)  #c(18, 16, 22, 13)  c(25,10,6, 28)
  colnames(data) <- c("up", "down")
  rownames(data) <- c("driver","suppressor")
  data;result
  # 执行费舍尔精确检验 
  fisher.test(result)  ;fisher.test(data)
  # 或者指定表格和检验方向的替代语法 fisher.test(x = data, alternative = "two.sided")
  #p-value = 0.1827
}

#output human significant ferroptosis related genes
df<-lisths(resulths,ferrdbhs,fc=0.25)
df

#8.0 immunemod免疫评分----
##8.1 install----
#https://github.com/cdesterke/immunemod
library(devtools);install_github("cdesterke/immunemod")
install_local("./data/immunemod-main")
##8.2 run----
library(immunemod)
data(data) #标准化exp_aver
data(phenotype)#分组信息
data(im)#基因集
res<-imenrich(exp_aver ,im,method="zscore",kcdf = "Gaussian")
head(res[,1:6])

phenotype<-data.frame(group=c(rep("Control",3),rep("Drug",3) ))
rownames(phenotype) <-colnames(exp_aver)

imheatmap(res,phenotype,scale="row",fontsize=10)
#图丑-自己画+显著性

library(pheatmap);library(RColorBrewer)
pheatmap(res[,1:6])
hallmark_p <-pheatmap(res,cluster_rows = T,cluster_cols  = F, show_colnames = T,
                      #scale="row",#标准化样本
                      border_color = "white",#cellwidth=10,cellheight = 10,
                      gaps_col= 3, #table(exp_marker1$response)
                      color= rev(colorRampPalette(brewer.pal(11, "Spectral"))(10) ) ,#colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
                      annotation_col=annotation_col#,annotation_row = anno_row
)
hallmark_p
save_pheatmap_pdf(hallmark_p, "./20240131out/pheatmap_immunemod-scalerow.pdf",5,3.5)

#显著性统计
pvalue_interest_res <-data.frame() #T.test
for (i in 1:nrow(res) ){
  print(i)
  pttest <- t.test(as.numeric(res[i,4:6] ),as.numeric(res[i,1:3]) )
  fc <- mean(as.numeric( res[i,4:6] ) )/mean( as.numeric(res[i,1:3]) )#drug/control
  pvalue_interest_res<-rbind(pvalue_interest_res,
                             data.frame(pvalue=pttest $p.value,
                                        FoldChange=fc,
                                        Symbol=rownames(res)[i] ) )
} 
#显著性差 inhibitory 0.03

#9.0 20250329- IB vs. DMSO wilcox----
#exp matrix
#romove duplicate
exp <-exp_aver

#样本分组
control <- colnames(exp)[grep("C", colnames(exp) )] #paired normal 
drug <- colnames(exp)[grep("D", colnames(exp) )] #tumor 



pvalue_g<-data.frame()

for (i in rownames(exp) ){ #[1:10]
  print(i)
  pwilcox <- t.test(as.numeric(exp[i,drug]),  #drug wilcox
                         as.numeric(exp[i,control]) ) #control
  fc <- mean( as.numeric(exp[i,drug]) )/mean( as.numeric(exp[i,control]) )#drug/control
  #fc <- as.numeric(pwilcox$estimate[1])/as.numeric(pwilcox$estimate[2])
  pvalue_g<-rbind(pvalue_g,data.frame(pvalue=pwilcox$p.value,FoldChange=fc,gene=i  ) )
  
}

pvalue_g$Sign. <- ifelse(pvalue_g$pvalue>= 0.05,"ns",
                         ifelse(pvalue_g$pvalue < 0.01,"**","*") )
pvalue_g$Fold <-ifelse(pvalue_g$FoldChange <=1,"< 1",
                       ifelse(pvalue_g$FoldChange >2,"> 2","> 1") )
rownames(pvalue_g) <-pvalue_g$gene

#write.csv(pvalue_g,"./out/IB_mRNA_Drug_Ctrl_wilcox.csv")
write.csv(pvalue_g,"./out/IB_mRNA_Drug_Ctrl_Ttest.csv")

