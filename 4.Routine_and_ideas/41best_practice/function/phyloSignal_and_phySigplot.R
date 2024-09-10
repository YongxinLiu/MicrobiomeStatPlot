
# #这个函数没有输出只会保存到指定的路径中，因为这个计算很浪费时间#---------
# phyloSignal(ps = ps,group  = "Group".path = "./")

phyloSignal = function(otu = NULL,
                       tax = NULL,
                       map = NULL,
                       tree = NULL ,
                       ps = NULL,
                       env = env,
                       group  = "Group",
                       path = "./"){

  # 抽平，默认使用最小序列抽平
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  # 抽平
  set.seed(72)
  psrare = rarefy_even_depth(ps)
  # 标准化
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))


  map = as.data.frame(sample_data(psrare))
  mapE =merge(map,env,by = "row.names",all= TRUE)
  row.names(mapE) = mapE$Row.names
  mapE$Row.names = NULL
  mapE$ID = row.names(mapE)
  sample_data(ps.norm) = mapE

  aa = levels(mapE$Group)

  dir.create(path)
  #----------分组计算门特尔相关,将结果保存，因为计算时间很长，只需计算一个就好了#-------
  eco = "Endosp."
  for (eco in as.character(unique(mapE$Group))){
    # Subset data
    print(paste("Now running", eco))
    # sub.physeq = phyloseq::subset_samples(ps.norm , Group == eco)
    sub.physeq = ps.norm
    otu = as.data.frame(vegan_otu(ps.norm))
    head(otu)
    map = as.data.frame(sample_data(ps.norm))
    mapsub <- map[map$Group == eco,]
    
    sample_data(sub.physeq) = mapsub
    
    # Remove OTUs not found in at least 3 samples
    OTU.table = otu_table(sub.physeq)
    OTU.table[OTU.table > 0] = 1
    OTU.freq = rowSums(OTU.table)
    OTU.freq = OTU.freq[OTU.freq > 2]
    sub.physeq = prune_taxa(names(OTU.freq), sub.physeq)
    sub.physeq

    # get phylogenetic distances
    tree = phy_tree(sub.physeq)
    phylo.dist = cophenetic(tree)
    sample_OTUs = tree$tip.label
    sam.phylo.dist = phylo.dist[sample_OTUs, sample_OTUs]
    sam.phylo.dist[upper.tri(sam.phylo.dist, diag=TRUE)] = NA

    # Generate dataframe of niche preference for pH, SOC and CN

    # site.chem.mat =  data.frame(sample_data(sub.physeq)) %>%
    #   # mutate(CN = percent_C / percent_N) %>%
    #   dplyr::select(ID, colnames(env))
    site.chem.mat =  env[row.names(env) %in% row.names(mapsub),]
    # rownames(site.chem.mat) = site.chem.mat$ID
    # site.chem.mat$ID = NULL
    site.chem.mat = as.matrix(site.chem.mat)

    otu.table = t(otu_table(sub.physeq))
    # head(otu.table)
    match(row.names(otu.table),row.names(site.chem.mat))
    
    OTU.niche = wascores(site.chem.mat, otu.table)
    OTU.niche.df = data.frame(OTU.niche)
    head( OTU.niche.df)
    # i =1
    for (i in 1:dim(OTU.niche.df)[2]) {
      pH.pref = OTU.niche.df[[i]]
      names(pH.pref) = rownames(OTU.niche.df)
      pH.dist = as.matrix(dist(pH.pref), labels=TRUE)
      sam.pH.dist = pH.dist[sample_OTUs, sample_OTUs]
      sam.pH.dist[upper.tri(sam.pH.dist, diag=TRUE)] = NA

      sam.pH.crlg = mantel.correlog(sam.pH.dist, sam.phylo.dist)
      # ?mantel.correlog
      filename = paste(path,eco,colnames(OTU.niche.df[i]), "_crlg.rds", sep="_")
      saveRDS(sam.pH.crlg, file=filename)

    }


  }

}






# #-----------提取数据#--------用于出图#-------
# #-本来我是不想分为两个函数的，但是由于上一步计算量太大，所以我分开
# #这里提供路径是为了读取上一步保存的文件
#
# result = phySigPlot(ps = ps,group  = "Group",path = "./")
#
# #提取图片
# p = result[[1]]
# p
# #-提取作图数据
# data = result[[2]]
# head(data)
# OTU.niche.df

# result = phySigPlot(ps = ps,group  = "Group",env = env,path = phypath)
#
# ps = ps
# group  = "Group"
# env = env
# path = "./"

phySigPlot = function(otu = NULL,
                      tax = NULL,
                      map = NULL,
                      tree = NULL,
                      ps = NULL,
                      group  = "Group",
                      env = env,
                      path = "./"){
  # 抽平，默认使用最小序列抽平
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  mapE = as.data.frame(sample_data(ps))
  for (eco in levels(mapE$Group)) {
      # eco = "KO"
      # i = 1
    for (i in 1:length(colnames(env))) {

      ag.pH.crlg = data.frame(readRDS(file=paste(path,eco,colnames(env[i]), "_crlg.rds", sep="_"))$mantel.res) %>%
        mutate(Group = eco, property = colnames(env)[i])

      if (i == 1) {

        data = ag.pH.crlg
      }
      if (i != 1) {

        data = rbind(data,ag.pH.crlg )
      }


    }

    if (eco == levels(mapE$Group)[1]) {
      data2 =  data
    }
    if (eco != levels(mapE$Group)[1]) {
      data2 = rbind(data2, data)
    }


  }

  dim(data2)




  eco.crlg = data2 %>%
    mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
    filter(!(is.na(Pr.corrected.)))
  eco.crlg$Group= factor(eco.crlg$Group)

  p = ggplot(data=eco.crlg, aes(x=class.index, y=Mantel.cor)) +
    geom_point(data=eco.crlg[eco.crlg$sig=="significant",], color = "black", size=2, shape=16) +
    geom_point(data=eco.crlg[eco.crlg$sig=="non-significant",], color = "black",size=2, shape=1) +
    geom_line(data=eco.crlg, aes(color=property)) +
    geom_hline(yintercept = 0, linetype=2) +
    labs(x = "Phylogenetic distance class", y="Mantel correlation", color="property") +
    # facet_grid(~Group)
    facet_wrap(~Group,scales="free_y",ncol  = 4)

  return(list(p,eco.crlg,data2))
}
