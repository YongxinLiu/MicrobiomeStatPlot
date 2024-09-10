




# #示例数据，详情加载 rdacca.hp 包后 ?mite
# library(rdacca.hp)
# library(vegan)
# library(ggClusterNet)
# library(tidyverse)
# 
# data(ps)
# data(env1)
# 
# mite = ps %>% 
#   filter_OTU_ps(200) %>%
#   vegan_otu()  %>% as.data.frame()
# mite.env =env1
# 
# hpPath = "./rdacca.hp/"
# dir.create(hpPath)
# 
# rdacca.hp.micro.p(
#   OTU = mite,
#   env = env1,
#   hpPath = hpPath,
#   cca = FALSE,
#   dbRDA = FALSE
# )




rdacca.hp.micro.p <- function(OTU = mite,
                              env = env1,
                              hpPath = hpPath,
                              cca  = FALSE,
                              dbRDA = FALSE,
                              rep = 9
){
  set.seed(123)
  permu_hp <- permu.hp(dv =OTU, iv = env, method = 'RDA', type = 'adjR2', permutations = rep)
  permu_hp
  #简单作图
  permu_hp$Variables <- rownames(permu_hp)
  permu_hp$p <- unlist(lapply(as.character(permu_hp$'Pr(>I)'), function(x) unlist(strsplit(x, ' '))[2]))
  #如需输出层次分割结果
  dat = permu_hp %>% as.data.frame() %>% arrange(desc(Individual))
  dat$id = row.names(dat)
  dat$id = factor(dat$id,levels = dat$id)
  
  p = ggplot(dat) + 
    geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) +
    geom_text(aes(x = id,y = Individual ,label = `Pr(>I)`), vjust = -0.3) +
    theme_classic()
  
  filename = paste(hpPath,'Mciro.rda.hp.p.csv',sep = "")
  write.csv(dat,filename)
  filename = paste(hpPath,'Micro.rda.hp.p.pdf',sep = "")
  ggsave(filename,p)
  filename = paste(hpPath,'Micro.rda.hp.p.jpg',sep = "")
  ggsave(filename,p)
  
  ###典范对应分析（CCA）
  
  #不同于 RDA，CCA 中的物种丰度一般不做 Hellinger 转化处理，可直接使用原始丰度矩阵
  #使用层次分割在 CCA 中分解每个环境变量的解释，详情 ?rdacca.hp
  #本示例计算校正后的 R2，在这里默认基于 1000 次置换获取对 CCA 中校正后 R2 的估计
  
  if (cca == TRUE) {
    set.seed(123)
    permu_hp <- permu.hp(dv =OTU, iv = env, method = 'CCA', type = 'adjR2', permutations = rep)
    
    #简单作图
    permu_hp$Variables <- rownames(permu_hp)
    permu_hp$p <- unlist(lapply(as.character(permu_hp$'Pr(>I)'), function(x) unlist(strsplit(x, ' '))[2]))
    #如需输出层次分割结果
    dat = permu_hp %>% as.data.frame() %>% arrange(desc(Individual))
    dat$id = row.names(dat)
    dat$id = factor(dat$id,levels = dat$id)
    
    p = ggplot(dat) + 
      geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) +
      geom_text(aes(x = id,y = Individual ,label = `Pr(>I)`), vjust = -0.3) +
      theme_classic()
    
    filename = paste(hpPath,'Mciro.cca.hp.p.csv',sep = "")
    write.csv(dat,filename)
    filename = paste(hpPath,'Micro.cca.hp.p.pdf',sep = "")
    ggsave(filename,p)
    filename = paste(hpPath,'Micro.cca.hp.p.jpg',sep = "")
    ggsave(filename,p)
  }
  
  
  ###基于距离的冗余分析（db-RDA），或称典范主坐标分析（CAP）
  
  # #db-RDA 需要计算群落相异指数，下文以 Bray-curtis 距离为例
  # #关于 β 多样性或相异指数的计算：https://mp.weixin.qq.com/s/Jwcz2zOwL7y2eu5U3zhAWQ
  # mite.bray <- vegdist(mite, method = 'bray')
  # #db-RDA 分析环境因子对甲螨物种丰度的影响
  # #细节可参考前文：https://mp.weixin.qq.com/s/KIhGjTL1Tzc-QL7Z03LT_g
  # mite.cap <- dbrda(mite.bray~., mite.env)
  # summary(mite.cap)  #db-RDA 概要
  # exp_adj <- RsquareAdj(mite.cap)$adj.r.squared * mite.cap$CCA$eig/sum(mite.cap$CCA$eig)  #获取校正后的 R2
  # cap1_exp <- paste('CAP1:', round(exp_adj[1]*100, 2), '%')
  # cap2_exp <- paste('CAP2:', round(exp_adj[2]*100, 2), '%')
  # plot(mite.cap, display = c('wa', 'cn'), type = 'n', xlab = cap1_exp, ylab = cap2_exp)  #db-RDA 的简单作图，只显示样本点和环境变量
  # text(mite.cap, display = 'cn', col = 'blue', cex = 0.8)
  # points(mite.cap, display = 'wa', pch = 19, cex = 1)
  
  #使用层次分割在 db-RDA 中分解每个环境变量的解释，详情 ?rdacca.hp
  #本示例计算校正后的 R2，db-RDA 默认使用 Ezekiel 公式计算调校正后的 R2
  if (dbRDA == TRUE) {
    set.seed(123)
    permu_hp <- permu.hp(dv =OTU, iv = env, method = 'dbRDA', type = 'adjR2', permutations = rep)
    
    #简单作图
    permu_hp$Variables <- rownames(permu_hp)
    permu_hp$p <- unlist(lapply(as.character(permu_hp$'Pr(>I)'), function(x) unlist(strsplit(x, ' '))[2]))
    #如需输出层次分割结果
    dat = permu_hp %>% as.data.frame() %>% arrange(desc(Individual))
    dat$id = row.names(dat)
    dat$id = factor(dat$id,levels = dat$id)
    
    p = ggplot(dat) + 
      geom_bar(aes(x = id,y = Individual ),stat = "identity",colour="black",fill="#9ACD32" ) +
      geom_text(aes(x = id,y = Individual ,label = `Pr(>I)`), vjust = -0.3) +
      theme_classic()
    filename = paste(hpPath,'Mciro.dbRDA.hp.p.csv',sep = "")
    write.csv(dat,filename)
    filename = paste(hpPath,'Micro.dbRDA.hp.p.pdf',sep = "")
    ggsave(filename,p)
    filename = paste(hpPath,'Micro.dbRDA.hp.p.jpg',sep = "")
    ggsave(filename,p)
  }
}


