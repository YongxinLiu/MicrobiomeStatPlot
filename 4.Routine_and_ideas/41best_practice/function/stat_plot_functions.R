  # 程序功能：R语言统计、绘图常用包和参数
  # Functions: Frequency used packages and parameters for statistics and plotting in R
  # 主要内容 Main steps: 
  # - 1. 常用包自动安装和加载
  # - 2. 绘图参数和函数
  # - 3. 统计函数



# 1. 常用包自动安装和加载

## 1.1 安装CRAN来源常用包

  site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
  # 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
  package_list = c("vegan", "reshape2", "ggplot2", "devtools", "bindrcpp", "ggthemes", "agricolae", "dplyr", 
                   "scales", "vegan", "pheatmap")
  # 判断R包加载是否成功来决定是否安装后再加载
  for(p in package_list){
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      install.packages(p, repos=site)
      suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
    }
  }

## 2.2 安装bioconductor常用包
  
  package_list = c("digest", "ggrepel")
  for(p in package_list){
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      source("https://bioconductor.org/biocLite.R")
      biocLite(p)
      suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
    }
  }

# 2.3 安装Github常用包

  package_list = c("kassambara/ggpubr")
  for(p in package_list){
    q=unlist(strsplit(p,split = "/"))[2]
    if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      install_github(p)
      suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
    }
  }


  
# 2. 绘图参数和函数
  
  # Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.

## 2.1 设置颜色
  
## 2.2 调置ggplot2主题
  
  
  alpha <- .7
  c_yellow <-          rgb(255 / 255, 255 / 255,   0 / 255, alpha)
  c_blue <-            rgb(  0 / 255, 000 / 255, 255 / 255, alpha)
  c_orange <-          rgb(255 / 255,  69 / 255,   0 / 255, alpha)
  c_green <-           rgb(  50/ 255, 220 / 255,  50 / 255, alpha)
  c_dark_green <-      rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
  c_very_dark_green <- rgb( 50 / 255, 150 / 255, 100 / 255, alpha)
  c_sea_green <-       rgb( 46 / 255, 129 / 255,  90 / 255, alpha)
  c_black <-           rgb(  0 / 255,   0 / 255,   0 / 255, alpha)
  c_grey <-            rgb(180 / 255, 180 / 255,  180 / 255, alpha)
  c_dark_brown <-      rgb(101 / 255,  67 / 255,  33 / 255, alpha)
  c_red <-             rgb(200 / 255,   0 / 255,   0 / 255, alpha)
  c_dark_red <-        rgb(255 / 255, 130 / 255,   0 / 255, alpha)
  
  # 主题为空白，轴线宽和颜色，刻度和文字，图例位置，正文文字
  main_theme = theme(panel.background=element_blank(),
                     panel.grid=element_blank(),
                     axis.line.x=element_line(size=.5, colour="black"),
                     axis.line.y=element_line(size=.5, colour="black"),
                     axis.ticks=element_line(color="black"),
                     axis.text=element_text(color="black", size=7),
                     legend.position="right",
                     legend.background=element_blank(),
                     legend.key=element_blank(),
                     legend.text= element_text(size=7),
                     text=element_text(family="sans", size=7))





#  3. 统计函数

# 函数summarySE：计算样本均值、标准差、标准误和置信区间
summarySE = function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                     conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # 计算长度
  length2 = function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac = ddply(data, groupvars, .drop=.drop,
                .fun = function(xx, col) {
                  c(N    = length2(xx[[col]], na.rm=na.rm),
                    mean = mean   (xx[[col]], na.rm=na.rm),
                    sd   = sd     (xx[[col]], na.rm=na.rm)
                  )
                },
                measurevar
  )
  # 重命名  
  datac = plyr::rename(datac, c("mean" = measurevar))
  # 计算标准偏差
  datac$se = datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult = qt(conf.interval/2 + .5, datac$N-1)
  datac$ci = datac$se * ciMult
  return(datac)
}


# vegan cca分析使用统计函数
variability_table <- function(cca){
  
  chi <- c(cca$tot.chi,
           cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  return(variability_table)
  
}

cap_var_props <- function(cca){
  
  eig_tot <- sum(cca$CCA$eig)
  var_propdf <- cca$CCA$eig/eig_tot
  return(var_propdf)
}

pca_var_props <- function(cca){
  
  eig_tot <- sum(cca$CA$eig)
  var_propdf <- cca$CA$eig/eig_tot
  return(var_propdf)
}

cca_ci <- function(cca, permutations=5000){
  
  var_tbl <- variability_table(cca)
  p <- permutest(cca, permutations=permutations)
  ci <- quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
  return(ci)
}


# 三元图绘制函数
tern_e<-function (x, scale = 1, dimnames = NULL, dimnames_position = c("corner",
                                                                       "edge", "none"), dimnames_color = "black", id = NULL, id_color = "black",
                  coordinates = FALSE, grid = TRUE, grid_color = "gray", labels = c("inside",
                                                                                    "outside", "none"), labels_color = "darkgray", border = "black",
                  bg = "white", pch = 19, cex = 1, prop_size = FALSE, col = "red",
                  main = "ternary plot", newpage = TRUE, pop = TRUE, ...)
{
  labels <- match.arg(labels)
  if (grid == TRUE)
    grid <- "dotted"
  if (coordinates)
    id <- paste("(", round(x[, 1] * scale, 1), ",", round(x[,
                                                            2] * scale, 1), ",", round(x[, 3] * scale, 1), ")",
                sep = "")
  dimnames_position <- match.arg(dimnames_position)
  if (is.null(dimnames) && dimnames_position != "none")
    dimnames <- colnames(x)
  if (is.logical(prop_size) && prop_size)
    prop_size <- 3
  if (ncol(x) != 3)
    stop("Need a matrix with 3 columns")
  if (any(x < 0))
    stop("X must be non-negative")
  s <- rowSums(x)
  if (any(s <= 0))
    stop("each row of X must have a positive sum")
  x <- x/s
  top <- sqrt(3)/2
  if (newpage)
    grid.newpage()
  xlim <- c(-0.03, 1.03)
  ylim <- c(-1, top)
  pushViewport(viewport(width = unit(1, "snpc")))
  if (!is.null(main))
    grid.text(main, y = 0.9, gp = gpar(fontsize = 18, fontstyle = 1))
  pushViewport(viewport(width = 0.8, height = 0.8, xscale = xlim,
                        yscale = ylim, name = "plot"))
  eps <- 0.01
  grid.polygon(c(0, 0.5, 1), c(0, top, 0), gp = gpar(fill = bg,
                                                     col = border), ...)
  if (dimnames_position == "corner") {
    grid.text(x = c(0, 1, 0.5), y = c(-0.02, -0.02, top +
                                        0.02), label = dimnames, gp = gpar(fontsize = 12))
  }
  if (dimnames_position == "edge") {
    shift <- eps * if (labels == "outside")
      8
    else 0
    grid.text(x = 0.25 - 2 * eps - shift, y = 0.5 * top +
                shift, label = dimnames[2], rot = 60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.75 + 3 * eps + shift, y = 0.5 * top +
                shift, label = dimnames[1], rot = -60, gp = gpar(col = dimnames_color))
    grid.text(x = 0.5, y = -0.02 - shift, label = dimnames[3],
              gp = gpar(col = dimnames_color))
  }
  if (is.character(grid))
    for (i in 1:4 * 0.2) {
      grid.lines(c(1 - i, (1 - i)/2), c(0, 1 - i) * top,
                 gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(1 - i, 1 - i + i/2), c(0, i) * top,
                 gp = gpar(lty = grid, col = grid_color))
      grid.lines(c(i/2, 1 - i + i/2), c(i, i) * top, gp = gpar(lty = grid,
                                                               col = grid_color))
      if (labels == "inside") {
        grid.text(x = (1 - i) * 3/4 - eps, y = (1 - i)/2 *
                    top, label = i * scale, gp = gpar(col = labels_color),
                  rot = 120)
        grid.text(x = 1 - i + i/4 + eps, y = i/2 * top -
                    eps, label = (1 - i) * scale, gp = gpar(col = labels_color),
                  rot = -120)
        grid.text(x = 0.5, y = i * top + eps, label = i *
                    scale, gp = gpar(col = labels_color))
      }
      if (labels == "outside") {
        grid.text(x = (1 - i)/2 - 6 * eps, y = (1 - i) *
                    top, label = (1 - i) * scale, gp = gpar(col = labels_color))
        grid.text(x = 1 - (1 - i)/2 + 3 * eps, y = (1 -
                                                      i) * top + 5 * eps, label = i * scale, rot = -120,
                  gp = gpar(col = labels_color))
        grid.text(x = i + eps, y = -0.05, label = (1 -
                                                     i) * scale, vjust = 1, rot = 120, gp = gpar(col = labels_color))
      }
    }
  xp <- x[, 2] + x[, 3]/2
  yp <- x[, 3] * top
  size = unit(if (prop_size)
    #emiel inserted this code. x are proportions per row.  x*s is original data matrix. s = rowsums of original data matrix (x*s)
    prop_size * rowSums(x*x*s) / max(  rowSums(x*x*s) )
    #prop_size * rowSums(    (x*s) * ((x*s)/s)) / max(  rowSums(    (x*s) * ((x*s)/s)) )
    else cex, "lines")
  grid.points(xp, yp, pch = pch, gp = gpar(col = col), default.units = "snpc",
              size = size, ...)
  if (!is.null(id))
    grid.text(x = xp, y = unit(yp - 0.015, "snpc") - 0.5 *
                size, label = as.character(id), gp = gpar(col = id_color,
                                                          cex = cex))
  if (pop)
    popViewport(2)
  else upViewport(2)
}


# 函数da_adonis：采用adonis对距离矩阵进行组间差异统计
# Compare each group distance matrix by vegan adonis in bray_curtis
da_adonis = function(sampleV){
  sampleA = as.matrix(sampleV$sampA)
  sampleB = as.matrix(sampleV$sampB)
  design2 = subset(sub_design, group %in% c(sampleA,sampleB))
  if (length(unique(design2$group))>1) {
    sub_dis_table = dis_table[rownames(design2),rownames(design2)]
    sub_dis_table = as.dist(sub_dis_table, diag = FALSE, upper = FALSE)
    adonis_table = adonis(sub_dis_table ~ group, data = design2, permutations = 1000) 
    adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
    print(paste("In", m, "pvalue between", sampleA, "and", sampleB, "is", adonis_pvalue, sep=" "))
    adonis_pvalue = paste(m, sampleA, sampleB, adonis_pvalue, sep="\t")
    return(adonis_pvalue)
  }
}


# 计算lm函数和评估R2，并格式化用于显示 
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
