rm(list=ls()) 
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


# 读取OTUtab或各级profiling.tab
Data = read.table("genus.profile.txt", header=T, row.names= 1, sep="\t", comment.char="") 
# 读取 design.txt or mapping file，第三列信息随意，没有回报错
design = read.table("design.txt", header=T, row.names= 1, sep="\t", comment.char="")

# 匹配design和Data的行名，用于进一步处理数据
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
# Acquisition mean profiling value of genus by groups (Sum also works)
Data4Pic = do.call(rbind, Data_mean)[-1,]
group = Data_mean$group
colnames(Data4Pic) = group
Data4Pic=as.data.frame(Data4Pic)
Data4Pic[Data4Pic>0]=1
# save data for further use
write.table(Data4Pic,"data4venn.txt",sep = "\t")

# Visualization with VennDiagram
# NOTICE: imagetype="pdf" will not work，only "npg","tiff",etc can be saved automatically by venn.plot
# BUT: we can save plot by define "filename=NULL", and save figure usibg "grid.draw"

pdf(file="Genus_venn.pdf", height = 4, width = 6)
p1 <- venn.diagram(
  x=list(
    A=row.names(Data4Pic[Data4Pic$A==1,]),#根据自己的分组，调整list中的分组情况
    B=row.names(Data4Pic[Data4Pic$B==1,]),
    C=row.names(Data4Pic[Data4Pic$C==1,]),
    D=row.names(Data4Pic[Data4Pic$D==1,])),
             filename = NULL, 
             lwd = 3,
             alpha = 0.6,
             label.col = "white",
             cex = 1.5,
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
             cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
             fontfamily = "serif",
             fontface = "bold",
             cat.fontfamily = "serif",
             cat.fontface = "bold",
             margin = 0.05)
p1
grid.draw(p1)
dev.off()

# Basic Visualization with UpsetR
pdf(file="Genus_Upsetplot.pdf", height = 4, width = 6)
p2 <-upset(Data4Pic, sets = colnames(Data4Pic),order.by = "freq")
p2
dev.off()

# Queries in UpsetR
pdf(file="Genus_Upsetplot_indiv.pdf",height = 4, width = 10)
p3<-upset(Data4Pic, sets = colnames(Data4Pic), mb.ratio = c(0.55, 0.45), order.by = "freq",
      queries = list(list(query=intersects, params=list("A", "B"), color="purple", active=T), 
                     list(query=intersects, params=list("C", "D", "A"), color="green", active=T), 
                     list(query=intersects, params=list("B", "C", "A", "D"), color="blue", active=T)), 
      nsets = 3, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Shared Genus",
      sets.x.label = "General Number in Each Group", text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5))
p3
dev.off()

# attribute.plots in UpsetR
library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
pdf(file="Genus_Upsetplot_attribute.pdf",height = 8, width = 10)
p4<-upset(Data4Pic, sets = colnames(Data4Pic), mb.ratio = c(0.55, 0.45), order.by = "freq",
          queries = list(list(query=intersects, params=list("A", "B"), color="purple", active=T), 
                         list(query=intersects, params=list("C", "D", "A"), color="green", active=T), 
                         list(query=intersects, params=list("B", "C", "A", "D"), color="blue", active=T)), 
          nsets = 3, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Shared Genus",
          sets.x.label = "General Number in Each Group", text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5), 
          attribute.plots = list(gridrows = 45, plots = list(list(plot = scatter_plot,  x = "ReleaseDate", y = "AvgRating", queries = T), 
                                                             list(plot = scatter_plot,  x = "AvgRating", y = "Watches", queries = F)), ncols = 2), 
          query.legend = "bottom")
p4
dev.off()
#### Upset参数解释 ####
# example：存放矩阵的变量名称；set：所需要的集合名称；mb.ratio：调整上下两部分的比例;
# order.by：排序方式，freq为按频率排序; queries：查询函数，用于对指定列添加颜色;
# param: list, query作用于哪个交集；color：每个query都是一个list，里面可以设置颜色,没设置的话将调用包里默认的调色板；
# active：active：被指定的条形图：TRUE显示颜色，FALSE在条形图顶端显示三角形;
# nset：集合数量，也可用set参数指定具体集合
# number.angles：上方条形图数字角度，0为横向，90为竖向，但90时不在正上方
# point.size：下方点阵中点的大小；line.size：下方点阵中每个线的粗细
# mainbar.y.label：上方条形图Y轴名称;sets.x.label：左下方条形图X轴名称
# text.scale：六个数字控制关系见；query.legend 指定query图例的位置

