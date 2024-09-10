
EdgerSuper2 = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,
                      ps = NULL,group  = "Group",pvalue = 0.05,
                      lfc =0,artGroup = NULL,
                      method = "TMM",
                      j = 2,
                      path = diffpath
){
  
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  # ps = ps %>% 
  #   ggClusterNet::tax_glom_wt(ranks = j)
  if (j %in% c("OTU","gene","meta")) {
    ps = ps 
  } else if (j %in% c(1:7)) {
    ps = ps %>% 
      ggClusterNet::tax_glom_wt(ranks = j)
  } else if (j %in% c("Kingdom","Phylum","Class","Order","Family","Genus","Species")){
    
  } else {
    ps = ps
    print("unknown j, checked please")
  }
  
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  colnames(sub_design) = "Group"
  Desep_group <- as.character(levels(as.factor(sub_design$Group)))
  Desep_group
  
  
  if ( is.null(artGroup)) {
    #--构造两两组合#-----
    aaa = combn(Desep_group,2)
    # sub_design <- as.data.frame(sample_data(ps))
  }
  if (!is.null(artGroup)) {
    aaa  = as.matrix(b )
  }
  otu_table = as.data.frame(ggClusterNet::vegan_otu(ps))
  count = as.matrix(otu_table)
  count <- t(count)
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  dim(sub_design)
  sub_design$SampleType = as.character(sub_design$Group)
  sub_design$SampleType <- as.factor(sub_design$Group)
  # create DGE list
  d = edgeR::DGEList(counts=count, group=sub_design$SampleType)
  d$samples
  d = edgeR::calcNormFactors(d,method=method)#默认为TMM标准化
  
  # Building experiment matrix
  design.mat = model.matrix(~ 0 + d$samples$group)
  colnames(design.mat)=levels(sub_design$SampleType)
  d2 = edgeR::estimateGLMCommonDisp(d, design.mat)
  d2 = edgeR::estimateGLMTagwiseDisp(d2, design.mat)
  fit = edgeR::glmFit(d2, design.mat)
  
  
  
  
  #------------根据分组提取需要的差异结果#------------
  for (i in 1:dim(aaa)[2]) {
    # i = 1
    Desep_group = aaa[,i]
    print( Desep_group)
    
    
    # head(design)
    # 设置比较组写在前面的分组为enrich表明第一个分组含量高
    
    group = paste(Desep_group[1],Desep_group[2],sep = "-")
    group
    BvsA <- limma::makeContrasts(contrasts =  group,levels=design.mat)#注意是以GF1为对照做的比较
    # 组间比较,统计Fold change, Pvalue
    lrt = edgeR::glmLRT(fit,contrast=BvsA)
    
    # FDR检验，控制假阳性率小于5%
    de_lrt = edgeR::decideTestsDGE(lrt, adjust.method="fdr", p.value=pvalue,lfc=lfc)#lfc=0这个是默认值
    summary(de_lrt)
    # 导出计算结果
    x=lrt$table
    x$sig=de_lrt
    head(x)
    #------差异结果符合otu表格的顺序
    row.names(count)[1:6]
    
    x <- cbind(x, padj = p.adjust(x$PValue, method = "fdr"))
    enriched = row.names(subset(x,sig==1))
    depleted = row.names(subset(x,sig==-1))
    
    x$level = as.factor(ifelse(as.vector(x$sig) ==1, "enriched",ifelse(as.vector(x$sig)==-1, "depleted","nosig")))
    x = data.frame(row.names = row.names(x),logFC = x$logFC,level = x$level,p = x$PValue)
    head(x)
    # colnames(x) = paste(group,colnames(x),sep = "")
    
    
    
    # x = res
    # head(x)
    #------差异结果符合otu表格的顺序
    # x = data.frame(row.names = row.names(x),logFC = x$log2FoldChange,level = x$level,p = x$pvalue) 
    x1 = x %>%
      dplyr::filter(level %in% c("enriched","depleted","nosig") )
    head(x1)
    x1$Genus = row.names(x1)
    # x$level = factor(x$level,levels = c("enriched","depleted","nosig"))
    if (nrow(x1)<= 1) {
      
    }
    x2 <- x1 %>% 
      dplyr::mutate(ord = logFC^2) %>%
      dplyr::filter(level != "nosig") %>%
      dplyr::arrange(desc(ord)) %>%
      head(n = 5)
    
    file = paste(path,"/",group,j,"_","Edger_Volcano_Top5.csv",sep = "")
    write.csv(x2,file,quote = F)
    head(x2)
    
    p <- ggplot(x1,aes(x =logFC ,y = -log2(p), colour=level)) +
      geom_point() +
      geom_hline(yintercept=-log10(0.2),
                 linetype=4,
                 color = 'black',
                 size = 0.5) +
      geom_vline(xintercept=c(-1,1),
                 linetype=3,
                 color = 'black',
                 size = 0.5) +
      ggrepel::geom_text_repel(data=x2, aes(x =logFC ,y = -log2(p), label=Genus), size=1) +
      scale_color_manual(values = c('blue2','red2', 'gray30')) + 
      ggtitle(group) + theme_bw()
    
    p
    
    file = paste(path,"/",group,j,"_","Edger_Volcano.pdf",sep = "")
    ggsave(file,p,width = 8,height = 6)
    
    file = paste(path,"/",group,j,"_","Edger_Volcano.png",sep = "")
    ggsave(file,p,width = 8,height = 6)
    
    
    colnames(x) = paste(colnames(x),sep = "")
    x$group = group
    
    
    if (i ==1) {
      table =x
    }
    if (i != 1) {
      table = rbind(table,x)
    }
  }
  x = table
  return(x)
  
}

