adaptive_lasso <- function(x, y, family = "binomial") {
  ridge3 <- glmnet(x = x, y = y, family = family, alpha = 0)
  ridge3_cv <- cv.glmnet(x = x, y = y, family = family, type.measure = "auc",  nfold = 5, alpha = 0)
  best_ridge_coef3 <- coef(ridge3_cv, s = ridge3_cv$lambda.min)
  best_ridge_coef3 <- as.numeric(best_ridge_coef3)[-1]
  
  alasso3 <- glmnet(x = x, y = y, family = family,
                    alpha = 1,
                    penalty.factor = 1 / abs(best_ridge_coef3))
  alasso3_cv <- cv.glmnet(x = x, y = y, family = family,
                          type.measure = "auc",
                          nfold = 5,
                          alpha = 1,
                          penalty.factor = 1 / abs(best_ridge_coef3),
                          keep = TRUE)
  
  m <- list(adalasso = alasso3, 
            adalasso_cv = alasso3_cv,
            adalasso_coef = coef(alasso3_cv, s = alasso3_cv$lambda.min))
  return(m)
}