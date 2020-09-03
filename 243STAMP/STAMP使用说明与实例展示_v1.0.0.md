# STAMP使用说明与实例展示



标签（空格分隔）： 微生物组 统计分析 可视化 STAMP


1. STAMP简介
----------

STAMP是一款分析微生物物种注释与功能注释的可视化软件，STAMP 1.0于2010年发表Bioinformatics，2014年的2.0版本同样在Bioinformatics发布，目前最新版本为2.1.3。截止到2020年8月15日，STAMP的引用次数分别达到了719次和1390次。该软件除了能够绘制探索性分析的相关图像之外，还提供了假设检验的功能。此外，STAMP采用了图形化界面，对用户比较友好。
![image.png-101.6kB][1]
![image.png-61.3kB][2]![image.png-66.6kB][3]

2. STAMP基本设置
------------

 - STAMP的输入文件

STAMP接受制表符分隔（tab-seprarated）的文件，也可以与主流生信软件如qiime、Mothur等对接（通过FIle--Create profile from... 实现）。文件包含注释层级和样本信息两部分，文件第一行为表头，含有注释信息的列应当是从最高级到最低级排列，且必须形成严格的树型结构。鉴于目前很多的分类分级系统（包括GreenGenes和SILVA等流行的分类法）的标签错误以及其他一些问题，STAMP网站提供了checkHierarchy.py这一脚本，可用于识别STAMP配置文件当中所有的非层级条目。而对于未知的条目，应记为uclassified（不区分大小写）。STAMP对于读取计数的形式没有特殊要求，可以为整数或任何实数，这使得标准化的方法可以不止一种。考虑到生物学数据低准确度、低精密度的特点，对于样本数量，STAMP的作者建议没有最小的样本数，具体的数量应当由样品本身决定，但如进行假设检验则必须符合相应的数据分布。
 ![image.png-35.7kB][4]
STAMP允许通过元数据（metadata）文件定义与样本相关联的其他数据。这一文件也应当是制表符分隔的文件。该文件的第一列表示每个样品的名称，并与STAMP配置文件当中的条目一一对应，其他列可以指定为与该样本相关的任何其他数据。
![image.png-20.7kB][5]

 - STAMP的假设检验

关于假设检验，STAMP提供了对多组、两组和两样品的统计检验方式，以及与之相应的事后检验 (Post-hoc test) 、置信区间宽度和多重检验等。对于多组、两组以及两样品的假设检验方法分别如下面表1、表2和表3所示。对于多组样品，作者推荐使用ANOVA进行假设检验，两组样品则建议使用Welch’s t-test这一适用性更广泛的检验方式，同时建议使用Fisher精确检验应对两样品比较的情况。多重检验校正方面，可以选择传统的Benjamini-Hochberg方法，但作者更偏向使用Storey’s FDR。这一方法的计算量更大，效果较Benjamini-Hochberg也更好。
![image.png-237.9kB][6]

表1：STAMP提供的对于多组样本的假设检验、事后检验与多重校正方法

![image.png-184.8kB][7]
![image.png-30.5kB][8]

表2：STAMP提供的对于两组样本的假设检验、置信区间与多重校正方法

![image.png-326.4kB][9]

表3：STAMP提供的对于两样品统计检验的情况所应用的假设检验、置信区间与多重检验校正方式

 

3. STAMP在微生物组研究当中的实例
--------------------

 - 例1：以不同频次使用抗菌药物的两组儿童唾液微生物的预测功能的变化

![image.png-225.5kB][10]

图4：在以低频次与高频次使用a) 全种类抗菌药物与b) 阿奇霉素的两组儿童的唾液微生物的MetaCyc功能预测。柱状图显示以PICRUST2预测的差异性MetaCyc通路的平均占比。组间差异显示95%的置信区间，并只显示Welch’s t-test经FDR校正后q value < 0.05的部分。

> Fig. 4 Functionally predicted MetaCyc pathways differing in proportions in high and low user groups of a) all AMs and in b) azithromycin. The bar plot shows mean proportions of differential MetaCyc pathways predicted using PICRUSt2. The difference in proportions between the groups is shown with 95% confidence intervals. Only p value < 0.05 (Welch’s t test, FDR adjusted), are shown and composition）

结果：在低频次和高频次使用抗菌药物的儿童当中，功能预测鉴定出21个显著差异的metaCyc通路（Fig. 4a）。这些通路在低频组中占比更高。差异最大的通路包括了L-精氨酸降解、L-谷氨酸降解、L-谷氨酸降解Ⅴ、多胺生物合成Ⅱ超通路以及嘌呤核苷酸降解Ⅱ。在使用阿奇霉素的两组当中，一共有十个差异通路 (Fig. 4b)。甲醇氧化至一氧化碳、L-精氨酸降解以及GDP-甘露糖生物合成通路在阿奇霉素的低频使用组中占比较高，而Kdo转移至脂质ⅣAⅢ、(5Z)-十二碳烯酸酯生物合成在高频使用组中占比更高。

> Functional predictions identified 21 differentially present metaCyc pathways between the low and high AM users when all AM use were combined (Fig. 4a). All of the pathways had higher proportions in the low AM use group. The largest significant differences were pathways for L-arginine degradation, L-glutamate degradation V, superpathway of polyamine biosynthesis II and purine nucleotides degradation II. Ten pathways differed between low and high azithromycin use (Fig. 4b). Methanol oxidation to carbon monoxide pathway, L-arginine degradation and GDP-mannose biosynthesis pathways showed higher proportions in the low azithromycin group, while Kdo transfer to lipid IVA III, (5Z)-dodecenoate biosynthesis and peptidoglycan maturation pathways showed higher proportions in the high azithromycin group

 - 例2：不同致病性大肠杆菌感染造成的微生物物种差异

![image.png-142.7kB][11]
FIG 5：黏附性弥散型大肠杆菌（DAEC）与肠毒性大肠杆菌（ETEC）感染中的丰度差异性物种。差异性物种的筛选条件为校正后p值小于等于0.05并且效应量（即组间差异大小）为0.8。（A和B）分别表示宏基因组分析当中注释为Fusobacterium mortiferum和Campylobacter concisus的序列所占百分比，(C 和 D)则分别为Bifidobacterium longum 和 Alloprevotella tannerae的。（E）为去除宿主与大肠杆菌的序列之后，根据宏基因组确定的分类组成（由MetaPhlAn2根据进化分支特异性标记基因注释到物种水平）所构建的PCA图。
 
> FIG 5 Differentially abundant (diagnostic) taxa between DAEC and ETEC infections. Differentially abundant species were reported if they had a corrected P value of ≤ 0.05 and an effect size (the magnitude of the difference between groups) of 0.8. (A and B) Proportions of metagenomic sequences assigned to Fusobacterium mortiferum and Campylobacter concisus, respectively. (C and D) Proportions of sequences assigned to Bifidobacterium longum and Alloprevotella tannerae, respectively. (E) PCA plot based on the taxonomic composition of each metagenome (annotated at the species level using clade-specific marker genes with MetaPhlAn2) after removal of human and E. coli reads from the libraries.

结果：对于DAEC和ETEC感染，在最初的物种注释当中至少有四个物种出现了差异。其中，Fusobacterium mortiferum (P = 0.025) 和Campylobacter concisus (P = 0.011) 在ETEC感染当中显著富集，而Bifidobacterium longum (P = 0.040) 和 Alloprevotella tannerae  (P = 0.046) 在DAEC感染当中丰度显著上升。基于物种水平的分类组成的PCA图显示ETEC感染的样品更相似，而DAEC组的样品则显示了更强的多样性。
> The initial taxonomic characterization revealed at least four species that were discriminatory of DAEC versus ETEC infections. Specifically, Fusobacterium mortiferum (P = 0.025) and Campylobacter concisus (P = 0.011) were significantly more abundant in ETEC infections (Fig. 5A and B), while Bifidobacterium longum (P = 0.040) and Alloprevotella tannerae (P = 0.046) were significantly more enriched in DAEC infections (Fig. 5C and D). A PCA based on taxonomic composition at the species level also revealed that metagenomes associated with ETEC infections tended to be taxonomically more similar among themselves, whereas DAEC samples showed more diversity.

 

4. STAMP实战
----------

 - 数据选择

这里选取STAMP使用手册当中的肠型数据，Enterotypes. profile. spf为分隔符分隔的OUT表，由门（Phyla）和属（Genera）两个分类层级构成；Enterotypes. metadata. tsv为tsv格式的元数据，由样本编号、肠型、国籍等信息组成。部分注释信息和样品元数据分别如下图所示。

![image.png-50kB][12]

![image.png-48.7kB][13]

 - 下载并安装软件

在浏览器地址栏输入https://beikolab.cs.dal.ca/software/STAMP ，在’Downloads’当中找到并点击STAMP v2.1.3 下载链接，保存安装程序。下载之后打开安装程序并选择路径进行安装，注意安装路径当中不得含有中文字符。
![image.png-152.8kB][14]
![image.png-22.6kB][15]
![image.png-65.1kB][16]

 - 多组比较

安装之后，打开STAMP，点击左上角的“file”-“load data”，分别导入Enterotypes. profile. spf和Enterotypes. metadata. tsv。
![image.png-19.1kB][17]
导入之后默认显示PCA结果，以散点图的形式展示门水平的差异：
![image.png-88.1kB][18]
点击“Configure plot”，设置图例位置于图像左上角，也可点击“View”-“Group legend” 查看分组信息。在‘Group field’选项当中重新分组，选择‘Enterotype’，并去除后三个非主要肠型，仅保留三种肠型。同时，更改“Profile level”为‘Genera’可以看到三种肠型在PCA图中分开较为明显。
![PCA configure plot.png-31.2kB][19]![PCA legend top left.png-43.3kB][20]
![PCA 3 enterotypes genera.png-103.7kB][21]


 - 切换图表类型

STAMP允许两组或多组样品以及两个样品之间的比较，支持的图标类型除了PCA图之外还有：

 a. 柱状图

显示每个样品当中feature的相对比例或序列数目（通过Configure plot设置），并添加组均值，图示为三种肠型当中拟杆菌属的相对丰度图。
![barplot.png-76.2kB][22]

 b. 箱线图

快速查看各组组内数据分布的基本情况，可通过'Show only active features'查看符合阈值的features。
![boxplot.png-10.4kB][23]

 c. 热图

显示每个Features在样品中丰度的比例，不仅显示所有样本的丰度值，还可以对行与列的各单元进行聚类显示之间的关系。通过选择'Show only active features'，可以看到三种肠型的样品有部分聚到一起，和PCA的结果较为接近。
![heatmap.png-49.7kB][24]

 d. Post-hoc 图

多组统计检验的无效假设(如ANOVA或Kruskal-Wallis)是所有组相等。提供每对组间测量的P-value和效应大小。在两组或两样品比较的情况下，Post-hoc 图则转换为Extended error bar，显示各feature在两组或两样品中的数据分布。
![post-hoc.png-11.6kB][25]
5. 参考资料
-------

STAMP. https://beikolab.cs.dal.ca/software/STAMP
STAMP User’s Guide. https://beikolab.cs.dal.ca/software/images/c/cd/STAMP_Users_Guide.zip
STAMP：扩增子、宏基因组统计分析神器(中文帮助文档). https://blog.csdn.net/woodcorpse/article/details/80458077
差异分析工具STAMP手册2:使用手册（汉版）. https://www.jianshu.com/p/331b6796f8ff
Parks DH, Tyson GW, Hugenholtz P, Beiko RG. STAMP: statistical analysis of taxonomic and functional profiles. Bioinformatics. 2014;30(21):3123-3124. doi:10.1093/bioinformatics/btu494
Parks DH, Beiko RG. Identifying biologically relevant differences between metagenomic communities. Bioinformatics. 2010;26(6):715-721. doi:10.1093/bioinformatics/btq041
Raju SC, Viljakainen H, Figueiredo RAO, et al. Antimicrobial drug use in the first decade of life influences saliva microbiota diversity and composition. Microbiome. 2020;8(1):121. Published 2020 Aug 21. doi:10.1186/s40168-020-00893-y
Peña-Gonzalez A, Soto-Girón MJ, Smith S, et al. Metagenomic Signatures of Gut Infections Caused by Different Escherichia coli Pathotypes. Appl Environ Microbiol. 2019;85(24):e01820-19. Published 2019 Nov 27. doi:10.1128/AEM.01820-19
Arumugam M, Raes J, Pelletier E, et al. Enterotypes of the human gut microbiome [published correction appears in Nature. 2011 Jun 30;474(7353):666] [published correction appears in Nature. 2014 Feb 27;506(7489):516]. Nature. 2011;473(7346):174-180. doi:10.1038/nature09944


  [1]: http://static.zybuluo.com/aldrich-cpu/kck8lljlwqua5ruedtannkdi/image.png
  [2]: http://static.zybuluo.com/aldrich-cpu/rl6xtume53hjz8eumgacvag6/image.png
  [3]: http://static.zybuluo.com/aldrich-cpu/j45atr31xl74zaolt7hz58z9/image.png
  [4]: http://static.zybuluo.com/aldrich-cpu/hmb3opmaar03ldmt99v8nym2/image.png
  [5]: http://static.zybuluo.com/aldrich-cpu/7p0ebvzgmauhf5elrchpfv9c/image.png
  [6]: http://static.zybuluo.com/aldrich-cpu/nbschbquu3qkboq8641dxd02/image.png
  [7]: http://static.zybuluo.com/aldrich-cpu/v5bnpjcyb5n7pbdna9ps067y/image.png
  [8]: http://static.zybuluo.com/aldrich-cpu/sc77p15lghlm1cs8s0cjw5wq/image.png
  [9]: http://static.zybuluo.com/aldrich-cpu/nok0pl7n5x8xm8pmlrz97vnn/image.png
  [10]: http://static.zybuluo.com/aldrich-cpu/89d26zkxd6q6fmwuxzz2l1zq/image.png
  [11]: http://static.zybuluo.com/aldrich-cpu/sisgld2u72zj0qvp105rez6g/image.png
  [12]: http://static.zybuluo.com/aldrich-cpu/3w77xvn2fw4nxyrusmfogpo7/image.png
  [13]: http://static.zybuluo.com/aldrich-cpu/89xfjlbilpjvvj3my1b3ee1u/image.png
  [14]: http://static.zybuluo.com/aldrich-cpu/y0bbq9cmxjhsev5i4wn7361d/image.png
  [15]: http://static.zybuluo.com/aldrich-cpu/y7rrnw1rv2uw9esxw3exew7k/image.png
  [16]: http://static.zybuluo.com/aldrich-cpu/8gyggcqhel370ryyhuv25a3a/image.png
  [17]: http://static.zybuluo.com/aldrich-cpu/po0v1jw80sgtirb6wj1ah069/image.png
  [18]: http://static.zybuluo.com/aldrich-cpu/ga6dlm8wrxm9ogloh110owju/image.png
  [19]: http://static.zybuluo.com/aldrich-cpu/2sq3re8q4nvmxhskr5zbcp9t/PCA%20configure%20plot.png
  [20]: http://static.zybuluo.com/aldrich-cpu/8losb51cpab7b6iui39wnp1y/PCA%20legend%20top%20left.png
  [21]: http://static.zybuluo.com/aldrich-cpu/a3lio51dxjnghctd0li6xrki/PCA%203%20enterotypes%20genera.png
  [22]: http://static.zybuluo.com/aldrich-cpu/5lbfagv3r0b0eouz69t601nt/barplot.png
  [23]: http://static.zybuluo.com/aldrich-cpu/72ioigwkrwfb5ten0tjzls8c/boxplot.png
  [24]: http://static.zybuluo.com/aldrich-cpu/5eewmgnhbvjj57xqz0r75384/heatmap.png
  [25]: http://static.zybuluo.com/aldrich-cpu/cs12thwsmd443jlsyte9l7ot/post-hoc.png