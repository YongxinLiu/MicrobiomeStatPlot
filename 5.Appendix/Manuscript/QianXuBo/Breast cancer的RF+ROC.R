library(mlbench)
library(caret)
library(ROCR)
data(BreastCancer)
bc_changed <- BreastCancer[2:11]  #仅选取一部分数据做分析

###把数据分为训练集和测试集
set.seed(100)
bc_rand <- bc_changed[order(runif(699)), ] #699 observations
bc_rand <- sample(1:699, 499) 
d_train <- bc_changed[ bc_rand,]
d_test  <- bc_changed[-bc_rand,]


###随机森林拟合
library(randomForest)
set.seed(100) 
rf <- randomForest(Class~., data=d_train, ntree=500, na.action=na.omit, importance=TRUE)

###ROC
prob <- predict(rf, type="prob", d_test)[,2]
mypred <- prediction(prob, d_test$Class)
myperformance <- performance(mypred, "tpr", "fpr")
plot(myperformance, col=2, colorize=T)

