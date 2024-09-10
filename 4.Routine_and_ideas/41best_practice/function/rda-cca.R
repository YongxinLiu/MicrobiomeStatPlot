

# otu = NULL
# tax = NULL
# map = NULL
# ps = ps
# env = envRDA
# path = RDApath
# group = "Group"
# chose.env = F


# library("Cairo")
# ggsave("geo_Fus_wilt.pdf", p1, width = 12, height =8 , device = cairo_pdf, family = "Song")

# ps = ps
# env = envRDA
# path = RDApath 

RDA_CCA = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,env = env,group = "Group",path = "./",chose.env = F){
  # library("vegan")
  # library("grid")
  dir.create(path)
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  
  # 相对丰度标准化编号：ps1_rela
  
  ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  mapping = as.data.frame( phyloseq::sample_data(ps_rela))
  env.dat = env
  samp.fg = colnames(otu)
  env.st = vegan::decostand(env.dat, method="standardize", MARGIN=2,na.rm = TRUE)#
  samp.env= rownames(env.st)
  my.env = match(samp.fg, samp.env)
  env.st3 = na.omit(env.st[my.env, ])  # omit the NA rows if without fg data
  samp.env= rownames(env.st3)
  my.fg = match(samp.env, samp.fg)
  otu = otu[, my.fg]
  ##without Latitude and Longitude
  env.st2=env.st3
  # for CCA calculation
  otu = t(otu)
  rowSums(otu)
  DCA= vegan::decorana(otu)  
  xxa = as.data.frame(DCA$rproj)
  if(max(xxa$DCA1)<=4 & max(xxa$DCA1)>=3){twochiose = "T"}else{twochiose = "F"}
  if(max(xxa$DCA1) > 4 | twochiose == "T") {
    ##choise CCA
    C.whole = vegan::cca(otu, env.st3)  ##rda(otu, env.st2)
    C.whole
    # for env selection by CCA inflation factors
    #Function vif.cca and alias.cca can be used to analyse linear dependencies among constraints and conditions.
    inf_factor = vegan::vif.cca(C.whole)
    
    # delete varable with max inflation factor
    na_env = which(is.na(inf_factor))
    if(isTRUE(length(na_env) > "0") ){
      inf_factor = inf_factor[-na_env]
    }
    
    max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
    env.st4 = env.st3
    while ( inf_factor[max_env] > 20){
      env.st4 = env.st4[,-max_env]
      C.reduced = vegan::cca(otu, env.st4)
      inf_factor = vegan::vif.cca(C.reduced)
      max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
    }
    output2 = inf_factor ;output2
    env.st4
    
    # for F and p values
    ind.p = array(0,dim=c(1,ncol(env.st4)))
    ind.F = array(0,dim=c(1,ncol(env.st4)))
    for(j in 1:ncol(env.st4)){
      ind.cca = vegan::cca(otu, env.st4[,j]) #ind.cca = cca(otu, env.st[,j], env.st[,-j])  #
      ind.sig = anova(ind.cca,step=1000)
      ind.p[1,j] = ind.sig$Pr[1]
      ind.F[1,j] = ind.sig$F[1]
    }
    
    colnames(ind.p) = colnames(env.st4)
    inf_Fp=rbind(output2,ind.F,ind.p)
    row.names(inf_Fp)=c("inf_factor","F","p")
    
    ##重新计算CCA
    C.whole = vegan::cca(otu, env.st4)  ##rda(otu, env.st3)
    x.sig = anova(C.whole)
    x.p = x.sig$Pr[1] ;x.p
    x.F = x.sig$F[1]  ;x.F
    F1 <- paste("anova F: ",round(x.F, 2), sep = "")
    pp1 = paste("p: ",round(x.p, 2), sep = "")
    title = paste(F1," ",pp1, sep = "")
    
    output1 = summary(C.whole)
    
    str(output1)
    a=output1$sites;a  ##样本坐标
    b=output1$cont$importance;b ##特征值，解释??? #eigenvals(C.whole)
    c=output1$biplot*5;c  ##环境因子坐标
    
    filenamea = paste(path,"cca_inf_Fp.txt",sep = "")
    write.table(inf_Fp,file=filenamea,sep="\t",col.names=NA)
    
    
    
    ca1=round(b[2,1],2);ca1
    ca2=round(b[2,2],2);ca2

    aa = merge(a,mapping, by="row.names",all=F);aa
    c = as.data.frame(c)
    
    library("ggplot2")
    library(ggrepel)
    # mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A","#E6AB02", "#B3DE69")
    p=ggplot()+
      geom_point(data=aa,aes(x=CCA1,y=CCA2, fill = Group),pch = 21,colour = "black",size = 4)+
      geom_segment(data = c,aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
                   arrow = arrow(angle=22.5,length = unit(0.2,"cm")),linetype=1, size=0.6,colour = "black")+
      # stat_ellipse( data=aa,linetype = 2,level = 0.65,aes(x=CCA1,y=CCA2,group  =SampleType, colour =  SampleType))+
      geom_text_repel(data = c,aes(x=CCA1,y=CCA2, label = row.names(c)))+
      
      labs(x=paste("CCA 1 (", ca1*100, "%)", sep=""),
           y=paste("CCA 2 (", ca2*100, "%)", sep=""),
           title=title)
    
    p
    p = p+theme_bw()+
      
      #scale_y_continuous(expand = c(0,0))+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
      theme(
        
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold")
        
        
      ) 
    p
    plotnamea = paste(path,"/CCA_env.pdf",sep = "")
    ggsave(plotnamea, p, width = 8, height = 6, device = cairo_pdf, family = "Song")
    
    plotnamea2 = paste(path,"/CCA_env.jpg",sep = "")
    ggsave(plotnamea2, p, width = 8, height = 6, device = cairo_pdf, family = "Song")
    
  } else if (max(xxa$DCA1) < 3| twochiose == "T" ){
    ##choise RDA
    C.whole = vegan::rda(otu, env.st2)
    C.whole
    
    

    if (chose.env == TRUE) {
      # for env selection by RDA inflation factors
      
      #Function vif.cca and alias.cca can be used to analyse linear dependencies among constraints and conditions.
      inf_factor = vegan::vif.cca(C.whole)
      inf_factor
      
      # delete varable with max inflation factor
      na_env = which(is.na(inf_factor))
      if(isTRUE(length(na_env) > "0") ){
        inf_factor = inf_factor[-na_env]
      }
      
      max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
      env.st4 = env.st3
      while ( inf_factor[max_env] > 20){
        env.st4 = env.st4[,-max_env]
        C.reduced = vegan::cca(otu, env.st4)
        inf_factor = vegan::vif.cca(C.reduced)
        max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
      }
      output2 = inf_factor ;output2
      
      
      
      ##重新计算rdA
      C.whole = vegan::rda(otu, env.st4)  ##rda(otu, env.st3)
    } else if (chose.env == FALSE) {
      
      C.whole = vegan::rda(otu, env.st2)
      C.whole
      inf_factor = vegan::vif.cca(C.whole)
      inf_factor
      output2 = inf_factor ;output2
      # for F and p values
      ind.p = array(0,dim=c(1,ncol(env.st2)))
      ind.F = array(0,dim=c(1,ncol(env.st2)))
      for(j in 1:ncol(env.st2)){
        ind.cca = vegan::cca(otu, env.st2[,j]) #ind.cca = cca(otu, env.st[,j], env.st[,-j])  #
        ind.sig = anova(ind.cca,step=1000)
        ind.p[1,j] = ind.sig$Pr[1]
        ind.F[1,j] = ind.sig$F[1]
      }
      
      colnames(ind.p) = colnames(env.st2)
      inf_Fp=rbind(output2,ind.F,ind.p)
      row.names(inf_Fp)=c("inf_factor","F","p")
    }
    
    
    
    x.sig = anova(C.whole)
    x.p = x.sig$Pr[1] ;x.p
    x.F = x.sig$F[1]  ;x.F
    F1 <- paste("anova F: ",round(x.F, 2), sep = "")
    pp1 = paste("p: ",round(x.p, 2), sep = "")
    title = paste(F1," ",pp1, sep = "")
    
    output1 = summary(C.whole)
    
    str(output1)
    a=output1$sites;a  ##样本坐标
    b=output1$cont$importance;b ##特征值，解释??? #eigenvals(C.whole)
    c=output1$biplot;c  ##环境因子坐标

    ca1=round(b[2,1],2);ca1
    ca2=round(b[2,2],2);ca2
    exp = c(ca1,ca2)
    aa = merge(a,mapping, by="row.names",all=F);aa
    c = as.data.frame(c)
    
    # library("ggplot2")
    # library(ggrepel)
    # mi=c("#1B9E77" ,"#D95F02", "#7570B3","#E7298A","#E6AB02", "#B3DE69")
    p=ggplot()+
      geom_point(data=aa,aes(x=RDA1,y=RDA2, fill = Group),pch = 21,colour = "black",size = 4)+
      geom_segment(data = c,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                   arrow = arrow(angle=22.5,length = unit(0.2,"cm")),linetype=1, size=0.6,colour = "black")+
      # stat_ellipse( data=aa,linetype = 2,level = 0.65,aes(x=RDA1,y=RDA2,group  =SampleType, colour =  SampleType))+
      ggrepel::geom_text_repel(data = c,aes(x=RDA1,y=RDA2, label = row.names(c)))+
      
      labs(x=paste("RDA 1 (", ca1*100, "%)", sep=""),
           y=paste("RDA 2 (", ca2*100, "%)", sep=""),
           title=title)
    
    p
    p = p+theme_bw()+
      
      #scale_y_continuous(expand = c(0,0))+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      # scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
      theme(
        
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold")
        
        
      ) 
    p

  }

  p3 = p+ ggrepel::geom_text_repel(data=aa, aes_string(x=colnames(aa)[2],y=colnames(aa)[3],label=colnames(aa)[1]),size=4)#?stat_el
  p3
  
  return(list(p,aa,p3,inf_Fp,c,title,exp))
  
}




RDA_CCA_explain_percent = function(ps = ps,env.dat = envRDA){
  #--变量相对丰度标准化编号：ps1_rela
  
  ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  mapping = as.data.frame( phyloseq::sample_data(ps_rela))
  
  
  # match env and fg datasets
  samp.fg = colnames(otu)
  env.st = vegan::decostand(env.dat, method="standardize", MARGIN=2,na.rm = TRUE)#
  samp.env= rownames(env.st)
  my.env = match(samp.fg, samp.env)
  env.st2 = na.omit(env.st[my.env, ])  # omit the NA rows if without fg data
  samp.env= rownames(env.st2)
  my.fg = match(samp.env, samp.fg)
  otu = otu[, my.fg]
  
  # for CCA calculation
  otu = t(otu)
  C.whole = vegan::rda(otu, env.st2)
  C.whole
  
  #----------------------------------------选择环境变量-----------------------------
  inf_factor = vegan::vif.cca(C.whole)
  inf_factor
  
  # delete varable with max inflation factor
  na_env = which(is.na(inf_factor))
  if(isTRUE(length(na_env) > "0") ){
    inf_factor = inf_factor[-na_env]
  }
  
  max_env = which(inf_factor == max(inf_factor,na.rm = T))
  env.st4 = env.st2
  
  while ( inf_factor[max_env] > 20){
    env.st4 = env.st4[,-max_env]
    C.reduced = vegan::cca(otu, env.st4)
    inf_factor = vegan::vif.cca(C.reduced)
    max_env = which(inf_factor == max(inf_factor,na.rm = T ))
  }
  
  
  output2 = inf_factor ;output2
  
  C.whole = vegan::rda(otu, env.st4)  ##rda(otu, env.st3)
  total.chi = C.whole$tot.chi
  ind.p = array(0,dim=c(1,ncol(env.st4)))
  for(j in 1:ncol(env.st4)){
    # j = 1
    ind.par = vegan::rda(otu, env.st4[,j], env.st4[,-j])
    ind.chi = ind.par$CCA$tot.chi
    ind.per = ind.chi/total.chi
    ind.p[j] = ind.per
  }
  ind.p
  
  rowname = colnames(env.st4);rowname
  out = matrix(data=NA,ncol=length(colnames(env.st4)),nrow=1);out
  out = ind.p
  rownames(out) = "percent"
  colnames(out) = rowname
  out
  
  
  
  #------提取解释比例
  total.chi = C.whole$tot.chi;total.chi
  total.constrained = C.whole$CCA$tot.chi ; total.constrained
  # 解释的比例
  explained.percent = (total.constrained) / total.chi;explained.percent
  # 未解释的比例
  unexplained.percent = (total.chi - total.constrained) / total.chi;unexplained.percent
  exp = data.frame(ID = c("explained.percent","unexplained.percent"),count = c(explained.percent,unexplained.percent))
  exp
  return(list(out,exp))
}

# result = RDA_CCA_explain_percent(ps = ps,env.dat = envRDA)
# 
# out = result[[1]]
# wxp = result[[2]]
# 
# filenamea = paste(RDApath,"each_env_exp_percent.csv",sep = "")
# write.csv(out,filenamea)
# 
# filenamea = paste(RDApath,"all_index_explain_percent.csv",sep = "")
# write.csv(exp,filenamea)
