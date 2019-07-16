mycv.glmnet <-
  function(x,
            y,
            standardize = FALSE,
            intercept = TRUE,
            penalty.factor = NULL,
            seed = NULL,
            nfolds = 5,
            family = "binomial",
            alpha = 1,
            measure = "auc",
            fun = "max",
            type = "response",
            weights = NULL,
            nlambda = 100,
            cores = 1,
            adaptive = FALSE) {
    
    predictors <- colnames(x)
    glm.formula <-  as.formula(paste('y ~', paste(colnames(x), collapse = '+')))
    
    if (family == 'binomial') {
      y <- as.factor(y)
    } 
    if (length(penalty.factor) == 0) {
      penalty.factor <- rep(1, ncol(x))
    } 
    if (is.null(weights)) {
      weights <- rep(1, nrow(x))
    } 
    if (!is.null(seed)) {
      set.seed(seed)
    }
    
    N <- length(y)

    foldsize <- round(N/ nfolds)
    fold.index <- list()
    fold.index[[nfolds]] <- NULL
    
    index <- c(1:N)
    
    for (i in 1:(nfolds-1)) {
      
      if (i == 1) {
        samp.index <- sample(index, size = foldsize, replace = FALSE)
        
      } else {
        samp.index <- sample(index[-unlist(fold.index)], size = foldsize, replace = FALSE)
      } 
      fold.index[[i]] <- samp.index
      
    } 
    fold.index[[nfolds]] <- index[unlist(fold.index)]
    
    glmnet.obj0 <- NULL
    
    if (adaptive == FALSE) {
      
      glmnet.obj0 <- glmnet(
        x,
        y,
        penalty.factor = penalty.factor,
        family = family,
        alpha = alpha,
        standardize = standardize,
        intercept = intercept,
        weights = weights,
        nlambda = nlambda
      )
      
      metrics0 <- matrix(rep(NA, nfolds*length(glmnet.obj0$lambda)), ncol = nfolds)
      
      for (i in 1:nfolds) {
        
        train.index <- unlist(fold.index[-i])
        test.index <- fold.index[[i]]
        
        train.glmnet <- glmnet(
          x = x[train.index, ],
          y  = y[train.index],
          penalty.factor = penalty.factor,
          family = family,
          lambda = glmnet.obj0$lambda,
          standardize = standardize,
          intercept = intercept,
          alpha = alpha,
          weights = weights[train.index],
          nlambda = nlambda
        )
        
        predict.y <- predict(train.glmnet, x[test.index,], type = type)
        
        metric.fold <- rep(NA, length(glmnet.obj0$lambda)) ## tu skonczylem
        
        for (j in 1:length(glmnet.obj0$lambda)) {
          if (family == "binomial") {
            metric.fold[j] <- glmnet::auc(y[test.index], 
                                          predict.y[, j], 
                                          weights[test.index])
            
          } else {
            metric.fold[j] <- sum(abs(predict.y[, j] - y[test.index]) * 
                                    (weights[test.index])) / sum(weights[test.index])
          } 
          } 
        metrics0[, i] <- metric.fold 
      } 
      
      metric0.mean <- apply(metrics0, 1, mean)
      
      best.metric0 <- eval(call(fun, metric0.mean))
      
      best.metric0.index <- which(metric0.mean == best.metric0)[1]
      
      return(
        list(
          coef = coef(glmnet.obj0)[, best.metric0.index],
          metrics = metrics0,
          lambda = glmnet.obj0$lambda,
          best.metric.index = best.metric0.index,
          method = 0,
          glmnet.obj0 = glmnet.obj0
        )
      ) 
    } 
    # adaptive 
    glmnet.obj1 <- NULL 
    glmnet.obj2 <- NULL 
    glmnet.obj3 <- NULL 
    glmnet.obj4 <- NULL 
    taus <- c(0.1,0.5,1,2); 
    tau1 <- taus[1]; 
    tau2 <- taus[2]; 
    tau3 <- taus[3]; 
    tau4 <- taus[4]; 
    # i n i t i a l weights, whole data
    
    data.ls.fit <-  data.frame(y=y, x); 
    ls.fit <-  NULL; 
    if (family == 'binomial' ) { 
      dsgn <-  svydesign(ids = ~1,
                         weights = weights,
                         data = data.ls.fit)
      
      ls.fit <- svyglm(glm.formula, design = dsgn, family = quasibinomial())
      
      }  else { 
      ls.fit <-  lm(glm.formula, data = data.ls.fit, weights = weights)
      
      }
    beta.init.whole <- (coef(ls.fit))[-1]; 
    
    # i n i t i a l weights, each f o l d 
    beta.list <- vector("list", nfolds)
    
    
    for (i in 1:nfolds) {
      train.index <- unlist(fold.index[-i])
      test.index <- fold.index[[i]];
      data.ls.fit  <- data.frame(y = y[train.index], x[train.index, ])
      ls.fit <- NULL

      if (family == "binomial") {
        dsgn <- svydesign(ids = ~1,
                          weights = weights[train.index],
                          data = data.ls.fit)
        ls.fit <- svyglm(glm.formula, design = dsgn, family = quasibinomial())
        
      } else {
        ls.fit <-
          lm(glm.formula, data = data.ls.fit , weights = weights[train.index])
        
      }
      beta.list[[i]] <- coef(ls.fit)[-1]; 
      }
      
    glmnet.objs <- list();
      
      for (i in 1:4) {
        glmnet.objs[[i]] <- glmnet(x, y, penalty.factor = (1/abs(beta.init.whole))^taus[i], 
                                   family = family , alpha = alpha , standardize = standardize ,
                                  intercept = intercept , weights = weights, nlambda = nlambda);
      }
      glmnet.obj1 <- glmnet.objs[[1]]; 
      glmnet.obj2 <- glmnet.objs[[2]]; 
      glmnet.obj3 <- glmnet.objs[[3]]; 
      glmnet.obj4 <- glmnet.objs[[4]];
      metrics <- vector("list", 4)
      
      for (k in 1:4) {
        metric.folds <- matrix(rep(NA, nfolds*length(glmnet.objs[[k]]$lambda)), ncol=nfolds );
        for (i in 1:nfolds) {
          train.index <- unlist(fold.index[-i]);
          test.index <- fold.index[[i]];
          beta.init <- beta.list[[i]];
          train.glmnet <- glmnet(x[train.index, ], 
                                y[train.index],
                                penalty.factor = (1/abs(beta.init))^taus[k], 
                                family = family, 
                                lambda = glmnet.objs[[k]]$lambda , 
                                standardize = standardize , intercept = intercept , alpha = alpha ,
                                weights = weights[train.index] , nlambda = nlambda);
          predict.y <- predict(train.glmnet, x[test.index, ], type = type); 
          metric.fold <- rep(NA, length(glmnet.objs[[k]]$lambda));
          grid.n <- min(length(train.glmnet$lambda), length(glmnet.objs[[k]]$lambda));
          for (j in 1:grid.n) {
            if (family == "binomial") {
              metric.fold[j] <- glmnet::auc(y[test.index], predict.y[, j], weights[test.index])
              
            } else {
              metric.fold[j] <- sum(abs(predict.y[, j] - y[test.index])*
                                      (weights[test.index])) / sum(weights[test.index])
              
            }
          } 
          metric.folds[, i] <- metric.fold 
          } 
        metrics[[k]] <- metric.folds
      } 


metrics1 <- metrics [[1]]
metrics2 <- metrics [[2]]
metrics3 <- metrics [[3]]
metrics4 <- metrics [[4]]
metric1.mean <- apply(metrics1, 1, mean, na.rm = T)
metric2.mean <- apply(metrics2, 1, mean, na.rm = T)
metric3.mean <- apply(metrics3, 1, mean, na.rm = T)
metric4.mean <- apply(metrics4, 1, mean, na.rm = T)
best.metric1 <- eval(call(fun, metric1.mean))
best.metric2 <-  eval(call(fun, metric2.mean))
best.metric3 <-  eval(call(fun, metric3.mean))
best.metric4 <-  eval(call(fun, metric4.mean))
best.metric <-  eval(call(fun, best.metric1, best.metric2, best.metric3, best.metric4))
best.metric.method <- which(c(best.metric1, best.metric2, best.metric3, best.metric4) == best.metric)[1]

if (best.metric.method == 1) {
  best.metric1.index <- which(metric1.mean == best.metric1)[1]
  return(
    list (
      coef = coef(glmnet.obj1)[, best.metric1.index],
      metrics = metrics1,
      lambda = glmnet.obj1$lambda,
      best.metric.index = best.metric1.index,
      method = 1,
      glmnet.obj = glmnet.obj1
    )
  )
  
} 
if (best.metric.method == 2) {
  best.metric2.index <- which(metric2.mean == best.metric2)[1] 
  return(
    list (
      coef = coef(glmnet.obj2)[, best.metric2.index],
      metrics = metrics2,
      lambda = glmnet.obj2$lambda,
      best.metric.index = best.metric2.index,
      method = 2,
      glmnet.obj = glmnet.obj2
      
    )
  )
} 
if (best.metric.method == 3) {
  best.metric3.index <- which(metric3.mean == best.metric3)[1] 
  return(
    list (
      coef = coef(glmnet.obj3)[, best.metric3.index],
      metrics = metrics3,
      lambda = glmnet.obj3$lambda,
      best.metric.index = best.metric3.index,
      method = 3,
      glmnet.obj = glmnet.obj3
    )
  )
} 
best.metric4.index = which(metric4.mean == best.metric4)[1] 
return(
  list (
    coef = coef(glmnet.obj4)[, best.metric4.index],
    metrics = metrics4,
    lambda = glmnet.obj4$lambda,
    best.metric.index = best.metric4.index,
    method = 4,
    glmnet.obj = glmnet.obj4
  )
)
}
