###属水平差异检验
library(ALDEx2)
library(ggplot2)
library(ggrepel)

#数据准备
mydata <- read.csv('feature_taxonomy_10.csv', skip=1, header=T)
mydata[1:5, 1:10]
genera <- mydata[,-c(1,2,3,4,6)]  #只剩下“属”
genera_sum <- aggregate(genera[ ,-1], by=list(genera[ ,1]), sum)
rownames(genera_sum) <- genera_sum[, 1]  #把“属”名作为行名
otu_g <- genera_sum[,-1]  #删除多余的第1列

#分组文件
group_letter_g <- substring(colnames(otu_g), 1, 1)
group <- ifelse(group_letter_g=="P", "JIA", "Control")

#统计过程
set.seed(1000)
otu_log_g <- aldex.clr(otu_g, group, mc.samples=128, verbose=TRUE)  #log转换，输入文件行必须是OTU，列是样本名
otu_test_g <- aldex.ttest(otu_log_g, group, paired.test=FALSE)  #Welch和Wilcoxon检验
head(otu_test_g)  #we.ep为welch的P值，we.eBH为校正P值；wi.ep为WilcoxonP值，wi.eBH为校正P值
otu_glm <- aldex.glm(otu_log_g, group)  #另一种检验方法

#Estimate Effect Size
otu_effect_g <- aldex.effect(otu_log_g, group, include.sample.summary=FALSE, verbose=FALSE)

#2个分析结果（P值文件和效应文件）合成一个文件
otu_g_all <- data.frame(otu_test_g, otu_glm[,3:4], otu_effect_g[,5:6])
head(otu_g_all)
#write.csv(otu_g_all, file="genera_p_values.csv")

#筛选出P值＜0.05的数据
subset(otu_g_all, otu_g_all$we.ep < 0.05)    #Welch检验P值＜0.05
subset(otu_g_all, otu_g_all$we.eBH < 0.05)   #Welch检验校正P值＜0.05
subset(otu_g_all, otu_g_all$wi.ep < 0.05)    #Wilcoxon检验P值＜0.05，P值未校正的情况下有8个属有差异
subset(otu_g_all, otu_g_all$wi.eBH < 0.05)   #Wilcoxon检验校正P值＜0.05，校正后有1个属有差异
subset(otu_g_all, otu_g_all$glm.ep < 0.05)   #glm检验P值＜0.05，P值未校正的情况下有13个属有差异
subset(otu_g_all, otu_g_all$glm.eBH < 0.05)  #glm检验校正P值＜0.05，校正后有2个属有差异：Dialister和Lachnospira
subset(otu_g_all, otu_g_all$effect < -0.5)   #查看有哪些属的effect＜-0.5





