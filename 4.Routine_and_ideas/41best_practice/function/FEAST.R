# source tracker analyses
#
# The function named 'FEAST'
#
#
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title FEAST Microbial traceability analysis
#' @description Input otutab, metadata or phyloseq object; ; output a table object.
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param group group ID;
#' @param sinkG object of sink group
#' @param sourceG object of source group
#' @details
#' By default, input phyloseq object include metadata and otutab
#' @return  a table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso beta_pcoa beta_cpcoa
#' @examples
#'
#' FEAST(otu = otu,map = map,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
#' FEAST(ps =ps,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
#'
#' @export
#'

# #清空内存
# rm(list=ls())
# #导入otu表格
# otu = read.delim("../../data/otutab.txt",row.names = 1)
# #导入分组文件
# map = read.delim("../../data/metadata.tsv",row.names = 1)
# head(map)
# result = FEAST(otu = otu,map = map,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
# result
# #-案例二
# ps = readRDS("../../data/ps_liu.rds")
# data(ps)
# result = FEAST(ps =ps,group = "Group",sinkG = "WT",sourceG = c("KO","OE"))
# result


FEAST = function(otu = otutab,map = metadata,ps = NULL,group = "Group",sinkG = "WT",sourceG = c("KO","OE"),
                 path = "E:/Shared_Folder/Function_local/R_function/micro/"
                 ){
  #-
  library(tidyverse)
  library("vegan")
  library("reshape2")
  # library(EasyMicrobiome)
  #

  source(paste(path,"FEAST-master/FEAST_src/src.R",sep = "/"))
  # #
  # if (is.null(ps) ) {
  #   head(otu)
  #   otu = as.matrix(otu)
  #   str(otu)
  #   colnames(map) = gsub(group,"AA", colnames(map))
  #   map$Group = map$AA
  #   map$Group = as.factor(map$Group )
  #
  #   ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE),
  #                  sample_data(map)
  #   )
  # }
  #
  # if (!is.null(ps) ) {
  #   ps = ps
  #   map = as.data.frame(sample_data(ps))
  #   map = map[, group]
  #   colnames(map) = "Group"
  #   map$Group = as.factor(map$Group)
  #   sample_data(ps) = map
  # }

  otu = NULL
  tax = NULL
  map = NULL
  Group = group
  ps = inputMicro(otu,tax,map,tree,ps,group  = Group)
  # 提取分组文件#----
  metadata <- as.data.frame(sample_data(ps))
  head(metadata)
  #--提取otu表格#--------
  otus <-  as.data.frame(t(vegan_otu(ps)))
  otus <- t(as.matrix(otus))


  head(metadata)
  #--将分组文件置于首列
  metadata$id = row.names(metadata)
  metadata = as.tibble(metadata)
  metadata<- dplyr::arrange(metadata, Group)


  #----提取样本名称，后续添加标记#------
  envs <- metadata$id

  #--目标ID提取#-----
  mu = metadata$id[metadata$Group==sinkG]


  # 设置FEAST运行参数#----
  EM_iterations = 1000 #default value
  different_sources_flag = 1
  Proportions_est <- list()

  #-提取每个分组测定的重复数量#-------
  rep = length(metadata$Group)/length(unique(metadata$Group))
  rep

  it = 1
  for(it in 1:rep){

    # it = 6
    #提取sink和source对应样本的位置，列的位置，方便后续提取#----
    train.ix <- which(metadata$Group%in%sourceG&metadata$id %in% metadata$id[seq(1, length(metadata$Group), rep)+(it-1)])
    test.ix <- which(metadata$Group==sinkG & metadata$id == mu[it])
    #--统计source样本数量#----
    num_sources <- length(train.ix)
    num_sources
    #-输入的是原始序列文件这里进行计算最小抽平数量#----------
    COVERAGE =  min(rowSums(otus[c(train.ix, test.ix),]))  #Can be adjusted by the user

    #提取source和sink对应的otu表格并抽平
    sources <- as.matrix(vegan::rrarefy(otus[train.ix,], COVERAGE))
    sinks <- as.matrix(vegan::rrarefy(t(as.matrix(otus[test.ix,])), COVERAGE))
    # sources = sources[1,1:10]

    if (length(sourceG) == 1) {
      tem1 <- sources %>% as.data.frame()
      sources = rbind(tem1,tem1) %>% as.matrix()
      dim(sources)
      row.names(sources) = paste(sourceG,1:2,sep = "")
    }

    #打印数量信息
    print(paste("Number of OTUs in the sink sample = ",length(which(sinks > 0))))
    print(paste("Seq depth in the sources and sink samples = ",COVERAGE))
    print(paste("The sink is:", envs[test.ix]))
    #FEAST 计算主函数#-------

    FEAST_output<-FEAST(source=sources,
                        sinks = sinks,
                        env = envs[train.ix],
                        em_itr = EM_iterations,
                        COVERAGE = COVERAGE)

    Proportions_est[[it]] <- FEAST_output$data_prop[,1]
    #整理结果#---
    names(Proportions_est[[it]]) <- c(as.character(envs[train.ix]), "unknown")

    if(length(Proportions_est[[it]]) < num_sources +1){
      tmp = Proportions_est[[it]]
      Proportions_est[[it]][num_sources] = NA
      Proportions_est[[it]][num_sources+1] = tmp[num_sources]
    }
    print("Source mixing proportions")
    print(Proportions_est[[it]])

  }

  went = as.data.frame(Proportions_est)
  colnames(went) = mu
  head(went)
  return(table = went)
}


#-------------
# source tracker result to plot
# The function named 'Plot_FEAST'
# You can learn more about package at:
#   https://github.com/microbiota/amplicon

#' @title FEAST Microbial traceability analysis result to plot
#' @description Input otutab, metadata or phyloseq object; ; output a table object.
#' @param result output object of the FEAST
#' @details
#' By default, input phyloseq object include metadata and otutab
#' @return  plot
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso beta_pcoa beta_cpcoa
#' @examples
#'
#'  Plot_FEAST(data = result)
#'
#'
#' @export
#'


# -例子
# Plot_FEAST(data = result)



Plot_FEAST = function(data = result){

  asx = as.data.frame(rowMeans(result))

  asx  = as.matrix(asx)
  asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100
  head(asx_norm)
  asx_norm = as.data.frame( asx_norm)
  colnames(asx_norm) = "present"
  # plotname = paste(path,"/FEAST_mean.pdf",sep = "")
  # pdf(file = plotname,width = 6,height = 6)
  labs <- paste0(row.names(asx_norm)," (", round(asx_norm[,1]/sum(asx_norm[,1])*100,2), "%)")
  asx_norm$lab = labs
  asx_norm$ID = row.names( asx_norm)
 p <-  ggplot(asx_norm, aes( x = "",y = present, fill = labs)) +
    geom_bar(stat = "identity",width = 20,color = "black") +
    coord_polar(theta = "y",direction=1) +
    theme_void()
 return(p)
}




# MuiPlot_FEAST(data = result)

# MuiPlot_FEAST = function(data = result){
#
#   par(mfrow=c(2,dim(result)[2]/2), mar=c(1,1,1,1))
#   # layouts = as.character(unique(metadata$SampleType))
#
#   for (i in 1:length(colnames(result))) {
#
#     labs <- paste0(row.names(result)," \n(", round(result[,i]/sum(result[,i])*100,2), "%)")
#
#     pie(result[,i],labels=labs, init.angle=90,col =  brewer.pal(nrow(result), "Reds"),
#         border="black",main =colnames(result)[i] )
#   }
# }
#
MuiPlot_FEAST = function(data = result){

  par(mfrow=c(2,dim(result)[2]/2), mar=c(1,1,1,1))
  # layouts = as.character(unique(metadata$SampleType))
  i = 1
  plots = list()
  for (i in 1:length(colnames(result))) {

    asx = data.frame(row.names = row.names(result),result[,i])

    asx  = as.matrix(asx)
    asx_norm = t(t(asx)/colSums(asx)) #* 100 # normalization to total 100
    head(asx_norm)
    asx_norm = as.data.frame( asx_norm)
    colnames(asx_norm) = colnames(result)[i]
    labs <- paste0(row.names(asx_norm)," (", round(asx_norm[,1]/sum(asx_norm[,1])*100,2), "%)")
    asx_norm$lab = labs
    asx_norm$ID = row.names( asx_norm)


    x <- colnames(result)[i]
    p <-  ggplot(asx_norm, aes( x = "",y = !!sym(x), fill = lab)) +
      geom_bar(stat = "identity",width = 20,color = "black") +
      coord_polar(theta = "y",direction=1) +
      labs(title = x) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(fill = guide_legend(title = NULL))
    p
    plots[[i]] = p
  }

  p  = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, legend="right")
  return(p)
}




