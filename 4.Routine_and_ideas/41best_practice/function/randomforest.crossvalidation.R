#For windows system#
##ramdomforest.crossvalidation.r##
##Begin##
rfcv1 <-
  function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5,
            mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE,
            ...)
  {
    classRF <- is.factor(trainy)
    n <- nrow(trainx)
    p <- ncol(trainx)
    if (scale == "log") {
      k <- floor(log(p, base = 1/step))
      n.var <- round(p * step^(0:(k - 1)))
      same <- diff(n.var) == 0
      if (any(same))
        n.var <- n.var[-which(same)]
      if (!1 %in% n.var)
        n.var <- c(n.var, 1)
    }
    else {
      n.var <- seq(from = p, to = 1, by = step)
    }
    k <- length(n.var)
    cv.pred <- vector(k, mode = "list")
    for (i in 1:k) cv.pred[[i]] <- rep(0,length(trainy))
    if (classRF) {
      f <- trainy
    }
    else {
      f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
    }
    nlvl <- table(f)
    idx <- numeric(n)
    for (i in 1:length(nlvl)) {
      idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold,
                                                  length = nlvl[i]))
    }
    res=list()
    for (i in 1:cv.fold) {
      all.rf <- randomForest(trainx[idx != i, , drop = FALSE],
                             trainy[idx != i],importance = TRUE)
      aa = predict(all.rf,trainx[idx == i, , drop = FALSE],type="prob")
      cv.pred[[1]][idx == i] <- as.numeric(aa[,2])
      impvar <- (1:p)[order(all.rf$importance[, 3], decreasing = TRUE)]
      res[[i]]=impvar
      for (j in 2:k) {
        12
        imp.idx <- impvar[1:n.var[j]]
        sub.rf <- randomForest(trainx[idx != i, imp.idx,
                                      drop = FALSE], trainy[idx != i]
        )
        bb <- predict(sub.rf,trainx[idx ==i,imp.idx, drop = FALSE],type="prob")
        cv.pred[[j]][idx == i] <- as.numeric(bb[,2])
        if (recursive) {
          impvar <- (1:length(imp.idx))[order(sub.rf$importance[,
                                                                3], decreasing = TRUE)]
        }
        NULL
      }
      NULL
    }
    if (classRF) {
      error.cv <- sapply(cv.pred, function(x) mean(factor(ifelse(x>0.5,1,0))!=trainy))
    }
    else {
      error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
    }
    names(error.cv) <- names(cv.pred) <- n.var
    list(n.var = n.var, error.cv = error.cv, predicted = cv.pred,res=res)
  }
##End##