[TOC]

## 前言
在微生物数据分析过程中，经常需要对某几组样本中共有或特有的OTU或微生物进行可视化。基于此需求，通常可以选择Venn图等进行可视化。然而，当分组信息过多，Venn的展示能力及可读性则有所下降，因此我们可以采用Venn图的升级版本——Upset Plot.

在本教程中，我将从菌群OTU表或者feature.tab开始，无需在代码外整理数据，从而直接实现用Venn或Upset Plot对集合信息进行可视化。进行可视化将使用VennDiagram和UpsetR等R包。

除了以上两个R包以外，还有ImageGP，yyplot,VennPainter，VennMaster以及TBtools的WonderfulVenn等软件或网站可供选择。详见：https://cloud.tencent.com/developer/article/1423035 

此外，需要指出的是，本文撰写过程中涉及到的描述和代码参考了诸多生信和可视化方面专家的前期工作，不具备开创性特色，不过搜集整理并让代码适用于直接对微生物分析中的OTU表格的分析。可以说本文的撰写，确实是站在巨人的肩上。文尾将对参考内容进行整理，此处不一一致谢了。

其他：本文所有示例数据和代码均已经上传Github：https://github.com/JerryHnuPKUPH/OTU2V-UpPlot

## Venn图简介

维恩图（英语：Venn diagram），又称Venn图。是十九世纪英国数学家约翰·维恩（John Venn）发明，用于展示集合之间大致关系的一类图形。其中圈或椭圆重合（overlap）的部分就是集合与集合间元素的交集，非重叠部分则为特定集合特有元素。

![image](https://bkimg.cdn.bcebos.com/pic/5d6034a85edf8db10a0e924f0b23dd54564e74fd?x-bce-process=image/resize,m_lfit,w_268,limit_1/format,f_jpg)


## Venn和Upset Plot对比

Venn和Upset Plot均可用于集合共有和特有元素信息进行可视化，但是当数据分组过多（>4），Venn图看起来会非常杂乱，而Upset plot可以展示≥5个以上分组的集合信息。

**总结起来：**

1）分组<5，Venn图更清晰；

2）分组≥5，Upset Plot更清晰；

3） Upset Plot展示方式更多元。


## 软件安装

```
# 设置清华源镜像
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list = c("VennDiagram","UpSetR")

# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
```
## VennDiagram和UpsetR的数据要求

值得注意的是，VennDiagram和UpsetR的数据要求并不相同，前者要求适用`list`输入以各组为集合的元素变量名,有几个分组就输入几个集合;而后者以元素变量名为行名，用数字0和1代表元素在分组集合中存在与否，数据输入是以数据框的形式输入。

**VennDiagram的数据类型**   --- | --- **UpsetR的数据类型**
| Set1 | Set2 | Set3 |      | ---  | Set1 | Set2 | Set3 |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| A    | A    | A    |      | A    | 1    | 1    | 1    |
| B    |      | B    |      | B    | 1    | 0    | 1    |
| C    | C    | C    |      | C    | 1    | 1    | 1    |
| D    | D    |      |      | D    | 1    | 1    | 0    |
| E    |      | E    |      | E    | 1    | 0    | 1    |
| F    |      |      |      | F    | 1    | 0    | 0    |



## 前期数据处理

```

# 读取OTUtab或各级profiling.tab
Data = read.table("genus.profile.txt", header=T, row.names= 1, sep="\t", comment.char="") 

# 读取 design.txt 和 mapping file，第三列信息随意，没有回报错
design = read.table("design.txt", header=T, row.names= 1, sep="\t", comment.char="")

# 匹配 design和 Data的行名，用于进一步处理数据
index = rownames(design) %in% colnames(Data) 
design = design[index,]
Data = Data[,rownames(design)] 
Data_t = t(Data)

# 合并design和Data数据表
Data_t2 = merge(design, Data_t, by="row.names")
# 删除来自design的非分组信息
Data_t2 = Data_t2[,c(-1,-3)]

# Define a function of mean
Data_mean = aggregate(Data_t2[,-1], by=Data_t2[1], FUN=mean)
# 计算各组数据的均值 (Sum值也行)
Data4Pic = do.call(rbind, Data_mean)[-1,]
group = Data_mean$group
# 为均值表添加列名（分组信息）
colnames(Data4Pic) = group
Data4Pic=as.data.frame(Data4Pic)

# 将数据表中>0的数值替换为1，数值>0则OTU或各级数据在分组中有出现，替换为1是用于可视化的数据准备
Data4Pic[Data4Pic>0]=1

# 保存数据，用于以后分析需要
write.table(Data4Pic,"data4venn.txt",sep = "\t")


```


## 使用VennDiagram绘制Venn图

**注意:** VennDiagram的图形绘制结果无法在Rstudio中直接呈现，而是直接生成图形并保存于工作目录。

**参数设置：**`x=list()`指定集合，由于VennDiagram要求输入以各组为集合的元素变量名，因此，我们将提取`Data4Pic`各组中数值`=1`的变量名作为数据输入的集合。
`filename=`指定图形绘制的结果保存的名称。
`imagetype=`参数设置图片生成的类型，但遗憾的是它只能指定`npg`,`tiff`等矢量图格式。
为了能够将图形绘制的结果保存为pdf格式，我们将`filename=`指定为`NULL`,并使用`grid.draw`函数输出图像。


```
pdf(file="Genus_venn.pdf") #设置pdf文件，用于存储绘制的图形
p1 <- venn.diagram(
  x=list(
    A=row.names(Data4Pic[Data4Pic$A==1,]),#提取各组中`=1`的行名。需根据自己的分组，调整list中的分组情况
    B=row.names(Data4Pic[Data4Pic$B==1,]),
    C=row.names(Data4Pic[Data4Pic$C==1,]),
    D=row.names(Data4Pic[Data4Pic$D==1,])),
             filename = NULL, #不指定保存的文件名，也不指定`imagetype`
             lwd = 3,
             alpha = 0.6,
             label.col = "white", #设置字体颜色
             cex = 1.5,
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"), #设置各组的颜色
             cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
             fontfamily = "serif", #设置字体
             fontface = "bold",
             cat.fontfamily = "serif",
             cat.fontface = "bold",
             margin = 0.05)
p1
grid.draw(p1) #使用`grid.draw`函数在`venn.diagram`绘图函数外绘制图形
dev.off()
```
图形绘制结果

![image](https://github.com/JerryHnuPKUPH/OTU2V-UpPlot/blob/master/Genus_venn.png?raw=true)

## 使用UpsetR绘制Upset Plot

```
# Upset plot的基本图形绘制
pdf(file="Genus_Upsetplot.pdf", height = 4, width = 6)
p2 <-upset(Data4Pic, sets = colnames(Data4Pic),order.by = "freq")
p2
dev.off()
```
**Upset参数解释1**
`example`：存放矩阵的变量名称；`set`：所需要的集合名称；`mb.ratio`：调整上下两部分的比例;`order.by`：排序方式，`freq`为按频率排序;

图形绘制结果

![image](https://github.com/JerryHnuPKUPH/OTU2V-UpPlot/blob/master/Genus_Upsetplot.png?raw=true)


```
# 使用Queries参数绘制Upset Plot
pdf(file="Genus_Upsetplot_indiv.pdf",height = 4, width = 10)
p3<-upset(Data4Pic, sets = colnames(Data4Pic), mb.ratio = c(0.55, 0.45), order.by = "freq",
      queries = list(list(query=intersects, params=list("A", "B"), color="purple", active=T), 
                     list(query=intersects, params=list("C", "D", "A"), color="green", active=T), 
                     list(query=intersects, params=list("B", "C", "A", "D"), color="blue", active=T)), 
      nsets = 3, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Shared Genus",
      sets.x.label = "General Number in Each Group", text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5))
p3
dev.off()

```
**Upset参数解释2**

`queries`：查询函数，用于对指定列添加颜色;`param: list`:query作用于哪个交集；`color`：每个query都是一个list，里面可以设置颜色,没设置的话将调用包里默认的调色板；`active`：被指定的条形图：`TRUE`显示颜色，`FALSE`在条形图顶端显示三角形;`nset`：集合数量，也可用set参数指定具体集合;`number.angles`：上方条形图数字角度，0为横向，90为竖向，但90时不在正上方;`point.size`：下方点阵中点的大小；`line.size`：下方点阵中每个线的粗细;`mainbar.y.label`：上方条形图Y轴名称;`sets.x.label`：左下方条形图X轴名称;`text.scale`：六个数字控制关系见；`query.legend`;指定query图例的位置…… 


图形绘制结果

![image](https://github.com/JerryHnuPKUPH/OTU2V-UpPlot/blob/master/Genus_Upsetplot_indiv.png?raw=true)

**更多参数信息**

此外，UpsetR还提供了`attribute.plots`参数,可绘制`histogram`,`scatter_plot`,`boxplot.summary`等图形的绘制。可以根据自己数据的需要进行配置数据，详细可见：https://cran.r-project.org/web/packages/UpSetR/vignettes/

**`attribute.plot`参数让UpsetR结果多元化**

*以下数据来源于`UpsetR`的示例数据，如需让自己的呈现的结果多元化，请结合自己的需求，进一步整理数据。*

（1）Upset Plot 与柱状图的组合

```
library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
    header = T, sep = ";")
upset(movies, main.bar.color = "black", queries = list(list(query = intersects, 
    params = list("Drama"), active = T)), attribute.plots = list(gridrows = 50, 
    plots = list(list(plot = histogram, x = "ReleaseDate", queries = F), list(plot = histogram, 
        x = "AvgRating", queries = T)), ncols = 2))
```

图形绘制结果

![image](https://i.postimg.cc/XY7TqVXs/7.png)

（2）Upset Plot与散点图的组合
```
upset(movies, main.bar.color = "black", queries = list(list(query = intersects, 
    params = list("Drama"), color = "red", active = F), list(query = intersects, 
    params = list("Action", "Drama"), active = T), list(query = intersects, 
    params = list("Drama", "Comedy", "Action"), color = "orange", active = T)), 
    attribute.plots = list(gridrows = 45, plots = list(list(plot = scatter_plot, 
        x = "ReleaseDate", y = "AvgRating", queries = T), list(plot = scatter_plot, 
        x = "AvgRating", y = "Watches", queries = F)), ncols = 2), query.legend = "bottom")
```
图形绘制结果

![image](https://i.postimg.cc/brgmft0p/8.png)

（3）Upset Plot与箱线图的组合

```
upset(movies, boxplot.summary = c("AvgRating", "ReleaseDate"))
```
图形绘制结果

![image](https://i.postimg.cc/bJb5WrFw/11.png)

(4) Upset Plot与自定义图的组合

```
# 自定义绘图参数`myplot`
myplot <- function(mydata, x, y) {
    plot <- (ggplot(data = mydata, aes_string(x = x, y = y, colour = "color")) + 
        geom_point() + scale_color_identity() + theme(plot.margin = unit(c(0, 
        0, 0, 0), "cm")))
}

another.plot <- function(data, x, y) {
    data$decades <- round_any(as.integer(unlist(data[y])), 10, ceiling)
    data <- data[which(data$decades >= 1970), ]
    myplot <- (ggplot(data, aes_string(x = x)) + geom_density(aes(fill = factor(decades)), 
        alpha = 0.4) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), legend.key.size = unit(0.4, 
        "cm")))
}

upset(movies, main.bar.color = "black", queries = list(list(query = intersects, 
    params = list("Drama"), color = "red", active = F), list(query = intersects, 
    params = list("Action", "Drama"), active = T), list(query = intersects, 
    params = list("Drama", "Comedy", "Action"), color = "orange", active = T)), 
    attribute.plots = list(gridrows = 45, plots = list(list(plot = myplot, x = "ReleaseDate", 
        y = "AvgRating", queries = T), list(plot = another.plot, x = "AvgRating", 
        y = "ReleaseDate", queries = F)), ncols = 2))
```
图形绘制结果

![image](https://i.postimg.cc/BQZhcms6/9.png)

## Reference

[1] https://www.cnblogs.com/jessepeng/p/11610055.html

[2] https://www.jianshu.com/p/285b4ac66768

[3] https://cran.r-project.org/web/packages/UpSetR/vignettes/

[4] https://blog.csdn.net/tuanzide5233/article/details/83109527