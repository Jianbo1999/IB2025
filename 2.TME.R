
#G:\IB\231-202306\4T1动物\DSP\测序结果\R\2.TME
if(T){
  rm(list = ls())
  setwd("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\R\\2.TME")
  #输出文件夹!!!!!!!
  program_name="TME_out"
  
  folder_path <- program_name#"./out" # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  #load exp and meta
  #exp<-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\1.Quality_Control\\raw_exp.csv",row.names = 1)
  #meta<-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\1.Quality_Control\\overall_annotation.csv",row.names = 1)
  #meta$id <-rownames(meta)#ID
  #meta$group <-meta$Group #比较分组信息!!!!!!!!!!
  load("G:/IB/231-202306/4T1动物/DSP/测序结果/R/2.TME/2.TME.RData")
}



#1.0 plot sample info----
df<-read.csv("sample_count_roi.csv")
library(ggplot2)
#柱状图
ggplot(df, aes(x = Type, y = Count,fill= Type)) + #,fill=Group
  geom_col() +
  scale_y_continuous(limits = c(0, max(df$Count)+3), expand = c(0, 0))+
  scale_fill_manual(values = c('#CCECFF',"skyblue","#4a6fe3") )+ #c("#66CCFF","#7EC7A7","#EDA065") c('#CCECFF',"skyblue","#4a6fe3")
  geom_text(aes(label = Count, vjust = -0.5),size = 3.5) +
  theme_bw()+#theme_minimal() +
  labs( x = "Area", y = "ROI count")+facet_grid(~Group)+
  theme(legend.position = "none",panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        axis.text = element_text(size = 12),axis.title = element_text(size = 12) )
ggsave(paste0(program_name,"/ROI_bar.pdf"),width = 2.5,height = 3,family="serif")


#百分柱状图
ggplot(df, aes(x = Group, y = Count,fill= Type)) + #,fill=Group
  geom_col(position = 'fill') +
  #scale_y_continuous(limits = c(0, max(df$Count)+3), expand = c(0, 0))+
  scale_fill_manual(values = c('#CCECFF',"skyblue","#4a6fe3") )+ #c("#66CCFF","#7EC7A7","#EDA065") c('#CCECFF',"skyblue","#4a6fe3")
  #geom_text(aes(label = Count, vjust = -0.5),size = 3.5) +
  theme_bw()+#theme_minimal() +
  labs( x = "Area", y = "Relative ROI count")+#facet_grid(~Group)+
  theme(panel.grid = element_blank(),
        #axis.text.x = element_text(angle = 60,hjust = 1),
        axis.text = element_text(size = 12),axis.title = element_text(size = 12) )
ggsave(paste0(program_name,"/ROI_bar_relative.pdf"),width = 2.5,height = 3,family="serif")
#手动添加count


#pie plot饼图
#https://www.jianshu.com/p/9963d9ad26df
library(ggplot2)

for (i in unique(df$Group) ){
  print(i)
  
  df1<-df[df$Group==i,]
  label <- paste0(df1$Count,'\n (', round(df1$Count/sum(df1$Count) * 100, 1), '%)')
  ggplot(df1, 
         aes(x = Group, y = Count,fill= Type)) +
    geom_bar(stat = 'identity', width = 0.5, position = 'stack')+ 
    scale_fill_brewer(palette ="Paired",direction = 1) +
    coord_polar(theta = 'y', direction = 1)+
    labs(x="",y="") +
    theme(axis.text = element_blank())+
    theme(axis.ticks = element_blank())+
    theme(panel.background = element_rect(I(0)))+
    geom_text(aes(x=1.1,label=as.character(label)),
              position = position_stack(reverse =F,vjust=0.5),size=3.5)
  ggsave(paste0(program_name,"/ROI_pie_",i,".pdf"),width = 3,height = 2.5,family="serif")
  
}

#2.0 TME:fraction-ImmCellAI----
#load result G:\IB\231-202306\4T1动物\DSP\测序结果\YKKY0396 交付\results\results\6.Signature\ImmCellAI
tme <-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\6.Signature\\ImmCellAI\\ImmCellAI_mouse_result.csv")
#合并meta
tme$id <-rownames(tme)
tme<-merge(meta,tme,by="id")#38:74
##2.1 DM-IB-all ----
###boxplot-all 合并柱状图----
library(reshape2)
bar_df <-melt(tme[c(38:73,34)],id.vars=c("Group") ) #exp_int, id.vars=c("Cell_type")
my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 

library(ggplot2);library(ggpubr)
p_method <-c("wilcox.test","t.test")
for (j in p_method){
  ggplot(bar_df,#
         aes(variable,value,color=Group ))+  #age_grou
    geom_point(alpha=0.5,size=0.2,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.7))+
    geom_boxplot(alpha=1,width=0.7,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    labs(x="",y="Cell proportion")+
    scale_color_manual(values =c("pink","skyblue") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")
    theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
    theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
          axis.text.x = element_text(angle = 90,hjust = 1),
          panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
    # stat_compare_means(#aes(group = Group) ,
    # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
    # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
    # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
    # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
    stat_compare_means(method = j,size=4, #"wilcox.test"
                       show.legend= F,label = "p.signif",#p.signif p.format
                       label.x =0.75,label.y =max(bar_df$value)-0.05)#
  ggsave( paste0(program_name,"/cibersort_boxplot_",j,".pdf"),width = 12,height = 6,family="serif")
  
}
###single-boxplot----
bar_df <-melt(tme[c(38:74,34)],id.vars=c("Group") ) #"Infiltration_score" 74
for (i in as.character(unique(bar_df$variable)) ){
  for (j in p_method[1] ){ #wilcox.test
    ggplot(bar_df[bar_df$variable %in% i,],#
           aes(Group,value,color=Group))+
      geom_boxplot(alpha=1,width=0.45,fill=NA,
                   position=position_dodge(width=0.5),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      geom_point(alpha=0.5,size=1,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.8))+
      # geom_violin(alpha=0.2,width=0.5,
      #             position=position_dodge(width=0.8),
      #             size=0.25)+
      labs(x="",y=i)+ #title = "",
      scale_color_manual(values =c("pink","skyblue","#9DC3E6") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
      theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
            #axis.text.x = element_text(angle = 90,hjust = 1),
            plot.title = element_text(size = 14),
            panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
      stat_compare_means(method = j,size=6, #"wilcox.test"
                         show.legend= F,label = "p.signif",#p.signif p.format
                         label.x =1.5,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20)
    #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
    ggsave( paste0(program_name,"/cibersort_boxplot_",i,"_",j,".pdf"),width = 2,height = 2.5,family="serif")#3 3.5
    
    print(i)
  }
}
###pheatmap----


##2.2 DM-IB-Type:TC IM TC_IM ----
for (x in unique(meta$Type) ){
  print(x)
  tme1 <- tme[tme$Type ==x,]
  
  ###Type 
  bar_df <-melt(tme1[c(38:73,34)],id.vars=c("Group") ) #exp_int, id.vars=c("Cell_type")
  my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  
  library(ggplot2);library(ggpubr)
  p_method <-c("wilcox.test","t.test")
  for (j in p_method){
    ggplot(bar_df,#
           aes(variable,value,color=Group ))+  #age_grou
      geom_point(alpha=0.5,size=0.2,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.7))+
      geom_boxplot(alpha=1,width=0.7,fill=NA,
                   position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      # geom_violin(alpha=0.2,width=0.9,
      #             position=position_dodge(width=0.8),
      #             size=0.25)+
      labs(x="",y="Cell proportion")+
      scale_color_manual(values =c("pink","skyblue") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")
      theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
            axis.text.x = element_text(angle = 90,hjust = 1),
            panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
      # stat_compare_means(#aes(group = Group) ,
      # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
      # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
      # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
      # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
      stat_compare_means(method = j,size=4, #"wilcox.test"
                         show.legend= F,label = "p.signif",#p.signif p.format
                         label.x =0.75,label.y =max(bar_df$value)-0.05)#
    ggsave( paste0(program_name,"/cibersort_boxplot_",j,"_",x,".pdf"),width = 12,height = 6,family="serif")
    
  }

  ###single boxplot
  bar_df <-melt(tme1[c(38:74,34)],id.vars=c("Group") ) #"Infiltration_score" 74
  for (i in as.character(unique(bar_df$variable)) ){
    for (j in p_method[1] ){ #wilcox.test
      ggplot(bar_df[bar_df$variable %in% i,],#
             aes(Group,value,color=Group))+
        geom_boxplot(alpha=1,width=0.45,fill=NA,
                     position=position_dodge(width=0.5),
                     size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                     outlier.stroke = 0.5)+
        geom_point(alpha=0.5,size=1,
                   position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                                 jitter.height = 0,
                                                 dodge.width = 0.8))+
        # geom_violin(alpha=0.2,width=0.5,
        #             position=position_dodge(width=0.8),
        #             size=0.25)+
        labs(x="",y=i)+ #title = "",
        scale_color_manual(values =c("pink","skyblue","#9DC3E6") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
        theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
        theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
              #axis.text.x = element_text(angle = 90,hjust = 1),
              plot.title = element_text(size = 14),
              panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
        stat_compare_means(method = j,size=6, #"wilcox.test"
                           show.legend= F,label = "p.signif",#p.signif p.format
                           label.x =1.5,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20)
      #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
      ggsave( paste0(program_name,"/cibersort_boxplot_",i,"_",j,"_",x,".pdf"),width = 2,height = 2.5,family="serif")#3 3.5
      
      print(i)
    }
  }
}
#not cibersort !just a name

#3.0 TME:immune function score----
#TC_IM IM 
#load geneset
geneset <-read.csv("./data/T_cell_function_geneset.csv")
#transform gene to mouse gene
library(biomaRt)
mg <- geneset$Genes#437 c("Spint1","Mgat4d","Wnt10a","Elmod2","Aspg")
#mouse to human

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#head(human@attributes[["name"]],20)
MtoH <- getLDS(attributes = "mgi_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "mgi_symbol", #参数过滤
               mart = mouse, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = mg, #要转换的基因集
               attributesL = c("hgnc_symbol"), #,"ensembl_gene_id" 要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = human, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE) #loss some genes ！！！！！！！！  mouse to human
#merge
colnames(MtoH)<-c("genes","Genes")
geneset_mouse <-merge(geneset,MtoH,by="Genes")
#与表达矩阵取交集
# common_gene <-intersect(geneset_mouse$genes,rownames(exp) )#178
# geneset_mouse<-geneset_mouse[geneset_mouse$genes %in% common_gene,]#310

geneset_mouse<-split(geneset_mouse$genes,geneset_mouse$Pathway)#geneset file input

geneset_mouse <-Filter(function(x) length(x) > 3, geneset_mouse )# >3个基因的基因集



#next time just use it in Rdata!

##3.1 TcellSI----
exp<-read.csv("G:\\IB\\231-202306\\4T1动物\\DSP\\测序结果\\YKKY0396 交付\\results\\results\\1.Quality_Control\\raw_exp.csv",row.names = 1)
#a function for TC_IM and IM type!

###
#devtools::install_github("GuoBioinfoLab/TCellSI")
#devtools::install_local("./data/TCellSI-main.zip")
library(TCellSI);library(ggplot2)
TcellSI <- CSS_Calculate(exp, 
                         ref = FALSE,#nbin = 1,
                         markers = geneset_mouse
                         )#genes

#Error in FUN(X[[i]], ...) : The number of features is less than 3
#未匹配到基因table(rownames(exp) %in% MtoH$genes ) #T 179
##
# Error in `ggplot2::cut_number()`:
#   ! Insufficient data values to produce 50 bins.
#过滤基因集>3无效；
#矩阵内容有问题？log2 exp无效
# boxplot(exp[1:100,])#0-15000
# boxplot(log2(exp[1:100,] ))#0-14
#0值太多？
library(dplyr)
exp1<- exp %>%
  filter(rowSums(.) > 95 ) #无法过滤，均高！#
#使用FPKM矩阵？？？？

#阳性对照成功！
boxplot(sample_expression[1:100,])#0-10
sample_expression <- TCellSI::exampleSample
TcellSI <- CSS_Calculate(sample_expression, 
                         ref = FALSE,#nbin = 10,
                         markers = split(geneset$Genes,geneset$Pathway)
                         )#genes


##3.2 ssGSEA----

#TC_IM IM boxplot
for (x in c("IM","TC_IM") ){
  print(x)
  meta1 <-meta[meta$Type ==x,]
  exp1 <-exp[,colnames(exp) %in% rownames(meta1) ]
  
  ####plot
  library(GSVA);library(GSEABase)
  gsva_matrix<- gsva(
    expr = as.matrix(exp1),
    gset.idx.list = geneset_mouse,
    method='ssgsea',
    kcdf='Gaussian'#,abs.ranking=TRUE
  ) #结果全部NA，使用全部基因表达
  gsva_matrix<-as.data.frame(t(gsva_matrix))
  gsva_matrix$id <-rownames(gsva_matrix)
  meta1 <-merge(meta1,gsva_matrix,by="id")
  ##single boxplot
  library(reshape2)
  bar_df <-melt(meta1[c(38:51,34)],id.vars=c("Group") ) #"Infiltration_score" 74
  for (i in as.character(unique(bar_df$variable)) ){
    for (j in p_method[1] ){ #wilcox.test
      library(ggplot2);library(ggpubr)
      ggplot(bar_df[bar_df$variable %in% i,],#
             aes(Group,value,color=Group))+
        geom_boxplot(alpha=1,width=0.45,fill=NA,
                     position=position_dodge(width=0.5),
                     size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                     outlier.stroke = 0.5)+
        geom_point(alpha=0.5,size=1,
                   position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                                 jitter.height = 0,
                                                 dodge.width = 0.8))+
        # geom_violin(alpha=0.2,width=0.5,
        #             position=position_dodge(width=0.8),
        #             size=0.25)+
        #labs(x="",y=i,title =x)+ #title = "",
        labs(x=x,y="ssGSEA score",title = i)+ #title = "",
        scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF"))+ # #"pink","skyblue","#9DC3E6") c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
        theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
        theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
              #axis.text.x = element_text(angle = 90,hjust = 1),
              plot.title = element_text(size = 14),
              panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
        stat_compare_means(method = j,size=4, #"wilcox.test"
                           show.legend= F,label = "p.format",#p.signif p.format
                           label.x =1.2,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20)
      #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
      #ggsave( paste0(program_name,"/ssgsea_boxplot_",i,"_",j,"_",x,".pdf"),width = 2,height = 2.5,family="serif")#3 3.5
      print(i)
    }
  }
  ##boxplot all
  for (j in p_method){
    ggplot(bar_df,#
           aes(variable,value,color=Group ))+  #age_grou
      geom_point(alpha=0.5,size=0.2,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.7))+
      geom_boxplot(alpha=1,width=0.7,fill=NA,
                   position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      # geom_violin(alpha=0.2,width=0.9,
      #             position=position_dodge(width=0.8),
      #             size=0.25)+
      labs(x="",y="ssGSEA score")+
      scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF") )+ #"pink","skyblue" #blue c("#0066CC","#66CCFF","#9DC3E6")
      theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
            axis.text.x = element_text(angle = 90,hjust = 1),
            panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
      # stat_compare_means(#aes(group = Group) ,
      # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
      # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
      # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
      # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
      stat_compare_means(method = j,size=4, #"wilcox.test"
                         show.legend= F,label = "p.signif",#p.signif p.format
                         label.x =0.75,label.y =max(bar_df$value)-0.05)#
    ggsave( paste0(program_name,"/ssgsea_boxplot_",j,"_",x,".pdf"),width = 6,height = 6,family="serif")
    
  }
  
  ####
}

#T_cell_co−inhibition 

##all sample
if(T){
  #exp meta #95
  
  ####plot
  library(GSVA);library(GSEABase)
  gsva_matrix<- gsva(
    expr = as.matrix(exp),
    gset.idx.list = geneset_mouse,
    method='ssgsea',
    kcdf='Gaussian'#,abs.ranking=TRUE
  ) #结果全部NA，使用全部基因表达
  gsva_matrix<-as.data.frame(t(gsva_matrix))
  gsva_matrix$id <-rownames(gsva_matrix)
  meta1 <-merge(meta,gsva_matrix,by="id")
  ##single boxplot
  library(reshape2)
  bar_df <-melt(meta1[c(38:51,34)],id.vars=c("Group") ) #"Infiltration_score" 74
  for (i in as.character(unique(bar_df$variable)) ){
    for (j in p_method[1] ){ #wilcox.test
      library(ggplot2);library(ggpubr)
      ggplot(bar_df[bar_df$variable %in% i,],#
             aes(Group,value,color=Group))+
        geom_boxplot(alpha=1,width=0.45,fill=NA,
                     position=position_dodge(width=0.5),
                     size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                     outlier.stroke = 0.5)+
        geom_point(alpha=0.5,size=1,
                   position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                                 jitter.height = 0,
                                                 dodge.width = 0.8))+
        # geom_violin(alpha=0.2,width=0.5,
        #             position=position_dodge(width=0.8),
        #             size=0.25)+
        #labs(x="",y=i,title =x)+ #title = "",
        labs(x="All",y="ssGSEA score",title = i)+ #title = "",
        scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF"))+ # #"pink","skyblue","#9DC3E6") c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
        theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
        theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
              #axis.text.x = element_text(angle = 90,hjust = 1),
              plot.title = element_text(size = 14),
              panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
        stat_compare_means(method = j,size=4, #"wilcox.test"
                           show.legend= F,label = "p.format",#p.signif p.format
                           label.x =1.2,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20)
      #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
      ggsave( paste0(program_name,"/ssgsea_boxplot_",i,"_",j,"_all.pdf"),width = 2,height = 2.5,family="serif")#3 3.5
      print(i)
    }
  }
  ##boxplot all
  for (j in p_method){
    ggplot(bar_df,#
           aes(variable,value,color=Group ))+  #age_grou
      geom_point(alpha=0.5,size=0.2,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.7))+
      geom_boxplot(alpha=1,width=0.7,fill=NA,
                   position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      # geom_violin(alpha=0.2,width=0.9,
      #             position=position_dodge(width=0.8),
      #             size=0.25)+
      labs(x="",y="ssGSEA score")+
      scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF") )+ #"pink","skyblue" #blue c("#0066CC","#66CCFF","#9DC3E6")
      theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
            axis.text.x = element_text(angle = 90,hjust = 1),
            panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
      # stat_compare_means(#aes(group = Group) ,
      # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
      # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
      # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
      # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
      stat_compare_means(method = j,size=4, #"wilcox.test"
                         show.legend= F,label = "p.signif",#p.signif p.format
                         label.x =0.75,label.y =max(bar_df$value)-0.05)#
    ggsave( paste0(program_name,"/ssgsea_boxplot_",j,"_all.pdf"),width = 6,height = 6,family="serif")
    
  }
  
  ####
}

#all is similar with IM TC_IM; T_cell_co−inhibition decrease !


#4.0 other immune inflitration ----
#https://zhuanlan.zhihu.com/p/709325517
# gmtFile=system.file("extdata", "SI_geneset.gmt",package="estimate")
# library(GSVA)
# library(limma)
# library(GSEABase)
# library(data.table)
# geneSet=getGmt(gmtFile,
#                geneIdType=SymbolIdentifier())
##4.1 CIBERSORT----
library(CIBERSORT)
Mice <- read.csv("./data/mice-signature.csv",row.names = 1)#https://www.nature.com/articles/srep40508#Sec23
Miceresults <- cibersort(sig_matrix = as.matrix(Mice), 
                         mixture_file = as.matrix(exp),
                         perm = 100,
                         QN = T)
Miceresults <-as.data.frame(Miceresults)
Miceresults$id <-rownames(Miceresults)
meta2<-merge(meta1,Miceresults[c(1:25,29)],by="id")#ssgesa+cibersort 52:76
##4.2 plot CIBERSORT----

#TC_IM IM boxplot
for (x in c("IM","TC_IM") ){
  print(x)
  meta1 <-meta2[meta$Type ==x,]
  #exp1 <-exp[,colnames(exp) %in% rownames(meta1) ]
  
  ##single boxplot
  library(reshape2)
  bar_df <-melt(meta1[c(52:76,34)],id.vars=c("Group") ) #"Infiltration_score" 74
  for (i in as.character(unique(bar_df$variable)) ){
    for (j in p_method[1] ){ #wilcox.test
      library(ggplot2);library(ggpubr)
      ggplot(bar_df[bar_df$variable %in% i,],#
             aes(Group,value,color=Group))+
        geom_boxplot(alpha=1,width=0.45,fill=NA,
                     position=position_dodge(width=0.5),
                     size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                     outlier.stroke = 0.5)+
        geom_point(alpha=0.5,size=1,
                   position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                                 jitter.height = 0,
                                                 dodge.width = 0.8))+
        # geom_violin(alpha=0.2,width=0.5,
        #             position=position_dodge(width=0.8),
        #             size=0.25)+
        #labs(x="",y=i,title =x)+ #title = "",
        labs(x=x,y="Cell fraction",title = i)+ #title = "",
        scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF"))+ # #"pink","skyblue","#9DC3E6") c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
        theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
        theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
              #axis.text.x = element_text(angle = 90,hjust = 1),
              plot.title = element_text(size = 14),
              panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
        stat_compare_means(method = j,size=4, #"wilcox.test"
                           show.legend= F,label = "p.format",#p.signif p.format
                           label.x =1.2,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20)
      #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
      #ggsave( paste0(program_name,"/CIBERSORT_boxplot_",i,"_",j,"_",x,".pdf"),width = 2,height = 2.5,family="serif")#3 3.5
      print(i)
    }
  }
  ##boxplot all
  for (j in p_method){
    ggplot(bar_df,#
           aes(variable,value,color=Group ))+  #age_grou
      geom_point(alpha=0.5,size=0.2,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.7))+
      geom_boxplot(alpha=1,width=0.7,fill=NA,
                   position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      # geom_violin(alpha=0.2,width=0.9,
      #             position=position_dodge(width=0.8),
      #             size=0.25)+
      labs(x="",y="Cell fraction")+
      scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF") )+ #"pink","skyblue" #blue c("#0066CC","#66CCFF","#9DC3E6")
      theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
            axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),
            panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
      # stat_compare_means(#aes(group = Group) ,
      # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
      # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
      # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
      # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
      stat_compare_means(method = j,size=4, #"wilcox.test"
                         show.legend= F,label = "p.signif",#p.signif p.format
                         label.x =0.75,label.y =max(bar_df$value)-0.05)#
    ggsave( paste0(program_name,"/REAL_CIBERSORT_boxplot_",j,"_",x,".pdf"),width = 8,height = 5,family="serif")
    
  }
  
  ####
}


##all sample
if(T){
  #exp meta #95
  
  ####plot

  ##single boxplot
  library(reshape2)
  bar_df <-melt(meta2[c(52:76,34)],id.vars=c("Group") ) #meta1 all sample types
  for (i in as.character(unique(bar_df$variable)) ){
    for (j in p_method[1] ){ #wilcox.test
      library(ggplot2);library(ggpubr)
      #bar_df[bar_df$variable %in% i,]$value <- scale(bar_df[bar_df$variable %in% i,]$value )
      bar <-bar_df[bar_df$variable %in% i,]
      bar$value <-scale(bar$value)
      
      ggplot(bar,#bar_df[bar_df$variable %in% i,]
             aes(Group,value, #value
                 color=Group))+
        geom_boxplot(alpha=1,width=0.45,fill=NA,
                     position=position_dodge(width=0.5),
                     size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                     outlier.stroke = 0.5)+
        geom_point(alpha=0.5,size=1,
                   position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                                 jitter.height = 0,
                                                 dodge.width = 0.8))+
        # geom_violin(alpha=0.2,width=0.5,
        #             position=position_dodge(width=0.8),
        #             size=0.25)+
        #labs(x="",y=i,title =x)+ #title = "",
        labs(x="All",y="Scaled cell fraction",#"Cell fraction"
             title = i)+ #title = "",
        scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF"))+ # #"pink","skyblue","#9DC3E6") c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
        theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
        theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
              #axis.text.x = element_text(angle = 90,hjust = 1),
              plot.title = element_text(size = 14),
              panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
        stat_compare_means(method = j,size=4, #"wilcox.test"
                           show.legend= F,label = "p.format",#p.signif p.format
                           label.x =1.2,label.y = max(bar$value)-max(bar$value)/20 #,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20
                           )
      #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
      #ggsave( paste0(program_name,"/CIBERSORT_boxplot_",i,"_",j,"_all.pdf"),width = 2,height = 2.5,family="serif")#3 3.5
      ggsave( paste0(program_name,"/CIBERSORT_boxplot_",i,"_",j,"_all_scaled.pdf"),width = 2,height = 2.5,family="serif")#3 3.5
      
      print(i)
    }
  }
  ##boxplot all
  for (j in p_method){
    ggplot(bar_df,#
           aes(variable,value,color=Group ))+  #age_grou
      geom_point(alpha=0.5,size=0.2,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.7))+
      geom_boxplot(alpha=1,width=0.7,fill=NA,
                   position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      # geom_violin(alpha=0.2,width=0.9,
      #             position=position_dodge(width=0.8),
      #             size=0.25)+
      labs(x="",y="Cell fraction")+
      scale_color_manual(values =c("#EDA065","#7EC7A7","#66CCFF") )+ #"pink","skyblue" #blue c("#0066CC","#66CCFF","#9DC3E6")
      theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
      theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
            axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),
            panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
      # stat_compare_means(#aes(group = Group) ,
      # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
      # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
      # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
      # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
      stat_compare_means(method = j,size=4, #"wilcox.test"
                         show.legend= F,label = "p.signif",#p.signif p.format
                         label.x =0.75,label.y =max(bar_df$value)-0.05)#
    #ggsave( paste0(program_name,"/REAL_CIBERSORT_boxplot_",j,"_all.pdf"),width = 8,height = 5,family="serif")
    
  }
  
  ####
}
##pheatmap----
#library(reshape2)
#bar_df <-melt(meta2[c(52:76,32:34)],id.vars=c("Group") ) #meta1 all sample types
#order
library(dplyr)
meta2 <-arrange(meta2,Group)
rownames(meta2)<-meta2$id 
#if numbers: Error in annotation_colors[[colnames(annotation)[i]]] :    subscript out of bounds

ann_colors <-list(Group=c("DM"="#EDA065","IB"="#7EC7A7"),
                  Type=c("IM" = "#A6CEE3", "TC" = "#CAB2D6","TC_IM" = "#B2DF8A" ) #,  #'Unchange'='#CAB2D6',
                  #Group=c(paste0(unique(group$group)[1]) ="#7EC7A7",paste(unique(group$group)[2]) ="skyblue") 
)
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

library(pheatmap);library(RColorBrewer)
hallmark_p<-pheatmap(meta2[,c(52:76)],scale = "row",#"column",#
         cluster_cols = T,cluster_rows = F,
         gaps_row = as.numeric(table(meta2$Group)[1]) ,
         color= colorRampPalette( c("navy", "white", "firebrick3"))(10),#rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)) ,#
         show_rownames  = F,border_color = NA,
         annotation_row = meta2[,c(32:34)],annotation_names_row =F,
         annotation_colors = ann_colors
)
hallmark_p
#save
save_pheatmap_pdf(hallmark_p, paste0(program_name,"/CIBERSORT_pheatmap.pdf"),5,5) #W H


hallmark_p<-pheatmap(t(meta2[,c(52:76)]),scale = "column",#"column",#row
         cluster_cols = F,cluster_rows = T,
         gaps_col = as.numeric(table(meta2$Group)[1]) ,
         color= colorRampPalette( c("navy", "white", "firebrick3"))(10),#rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)) ,#
         show_colnames  = F,border_color = NA,
         annotation_col = meta2[,c(32,34)],annotation_names_col =F,
         annotation_colors = ann_colors
         )
hallmark_p
save_pheatmap_pdf(hallmark_p, paste0(program_name,"/CIBERSORT_pheatmap-1.pdf"),5,4) #W H
rm(hallmark_p)

##add sign.添加显著性
score_df<-meta2[,c(52:76)] #as.data.frame(hallmark50_score)
#score_df<-scale(score_df)
sample_control<-rownames(meta2[meta2$Group=="DM",])
sample_drug<-rownames(meta2[meta2$Group=="IB",])

pvalue_g<-data.frame() #T.test?
for (i in colnames(score_df) ){
  print(i)
  pwilcox <- wilcox.test(as.numeric(score_df[sample_drug,i]),  #IB
                    as.numeric(score_df[sample_control,i]) ) #DM
  fc <- mean( as.numeric(score_df[sample_drug,i]) )/mean( as.numeric(score_df[sample_control,i]) )#drug/control
  #fc <- as.numeric(pwilcox$estimate[1])/as.numeric(pwilcox$estimate[2])
  pvalue_g<-rbind(pvalue_g,data.frame(pvalue=pwilcox$p.value,FoldChange=fc,gene=i  ) )
} #
#scale后为负值，且无差异

pvalue_g$Sign. <- ifelse(pvalue_g$pvalue>= 0.05,"ns",
                         ifelse(pvalue_g$pvalue < 0.01,"**","*") )
pvalue_g$FC <-ifelse(pvalue_g$FoldChange <=1,"< 1",
                     ifelse(pvalue_g$FoldChange >2,"> 2","> 1") )
rownames(pvalue_g) <-pvalue_g$gene

annotation_row_p <- data.frame(
  Sign.=pvalue_g$Sign.,FoldChange=pvalue_g$FC) #, size = 20, replace = TRUE  adj.P.Val
rownames(annotation_row_p) <- rownames(pvalue_g)
ann_colors_p=list(FoldChange=c('< 1'='#CAB2D6','> 1'="#FDBF6F",'> 2'="#B15928"),
                  Sign.=c("ns"="grey","*"="skyblue",'**'="blue" ),#""*"="#B2DF8A",'**'="#33A02C"
                  Group=c("DM"="#EDA065","IB"="#7EC7A7"),
                  Type=c("IM" = "#A6CEE3", "TC" = "#CAB2D6","TC_IM" = "#B2DF8A"))

library("gplots");library(viridis)

hallmark_g <-pheatmap(t(score_df),scale = "column",#"column",#row
         cluster_cols = F,cluster_rows = T,
         gaps_col = as.numeric(table(meta2$Group)[1]) ,
         color= colorRampPalette( c("navy", "white", "firebrick3"))(10),#rev(colorRampPalette(brewer.pal(11, "Spectral"))(8)) ,#
         show_colnames  = F,border_color = NA,
         annotation_col = meta2[,c(32,34)],annotation_names_col =F,
         annotation_row = annotation_row_p,annotation_names_row = F,
         annotation_colors = ann_colors_p#ann_colors
)
hallmark_g
save_pheatmap_pdf(hallmark_g, paste0(program_name,"/CIBERSORT_pheatmap-sign.pdf"),5.5,3.5) #W H
rm(hallmark_g)
