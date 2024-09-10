# Beta diversity statistics in all groups
#
# The function named 'MicroTest'
# which call adonis, anosim or MRPP to test beta-diversity. It usually called by BetaDiv.
#
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title Beta diversity statistics by adonis/anosim/MRPP in all groups
#' @description Input phyloseq object, test method and distance type
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param ps alternative input;
#' @param tree tree/nwk file;
#' @param dist distance type, including "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski"  "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial"  "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co";
#' @param group group ID;
#' @param method DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA;
#' @param pvalue.cutoff Pvalue threshold, default in 0.05;
#' @param Micromet statistics default by adonis, alternative anosim or MRPP;

#' @details
#' By default, input phyloseq object include metadata and otutab
#' The available diversity indices include the following:
#' \itemize{
#' \item{most used indices: bray unifrac wunifrac}
#' \item{other used indices: dpcoa jsd manhattan euclidean canberra kulczynski jaccard gower altGower morisita horn mountford raup binomial chao cao w -1 c wb r I e t me j sor m -2 co}
#' }
#' @return stat table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @seealso BetaDiv beta_pcoa beta_cpcoa
#' @examples
#' # Input otutab, metadata, and options group, method and distance
#' MicroTest(otu = otutab_rare, map = metadata, group = "Group", Micromet = "MRPP", dist = "bray")
#' # Input phyloseq format input, and options group, method and distance
#' MicroTest(ps = ps, group = "Group", Micromet = "MRPP", dist = "bray")
#' @export


MicroTest = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,
                     group = "Group", Micromet = "MRPP", dist = "bray"){

    
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  

  # 求取相对丰度
  
  ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela

  # library(vegan)
  # -准备矩阵和分组文件
  map = as.data.frame(phyloseq::sample_data(ps1_rela))

  unif = phyloseq::distance(ps, method=dist)

# adonis#----
  if (Micromet == "adonis") {
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    
    # a = round(as.data.frame(ado$aov.tab[5])[1,1],3)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    # b = as.data.frame(ado$aov.tab[6])[1,1]
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }

# MRPP#----
  if (Micromet == "MRPP") {
    dat.mrpp = vegan::mrpp(unif, map$Group)
    a = round(dat.mrpp$delta,3)
    R2 = paste("MRPP.delta ",a, sep = "")
    p_v = paste("p: ",round(dat.mrpp$Pvalue,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }

# anosim#----
  if (Micromet == "anosim") {
    dat.ano = vegan::anosim(unif, map$Group)
    a = round(dat.ano$statistic,3)
    R2 = paste("ANOSIM.r ",a, sep = "")
    p_v = paste("p: ",round(dat.ano$signif,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }
  return(title1)
}
