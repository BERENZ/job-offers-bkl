library(glmnet);
library(ROCR);
library(boot);
library(sampling);
library(survey);
library(caret);

options(warn=-1);

# Bootstrap with replacement
T.y2.lasso <- function(data, indices, sim.name, alpha, nlambda, nfolds, adaptive, method) {
  data.samp <- data[indices,];
  wts <- 1/data.samp$pk;
  y2.bin.samp <- data.samp[,"y2"];

  # LASSO with 
  X.lasso.y2.samp <- model.matrix(model,data=data.samp);

  That.y2.bin=NA;
  possibleError=0;
  tryCatch({

  cvobj.lasso.y2.bin <- mycv.glmnet(X.lasso.y2.samp[,-1],y2.bin.samp,family="binomial",seed=sum(indices),nfolds=nfolds,alpha=alpha,measure="auc",type="response",fun="max",weights=wts,nlambda=nlambda);

  coef.lasso.y2.bin <- cvobj.lasso.y2.bin$coef;
  lasso.fit.pop <- expit(cbind(X.y2.bench.lasso) %*% as.vector(coef.lasso.y2.bin));
  lasso.fit.samp <- expit(cbind(X.lasso.y2.samp) %*% as.vector(coef.lasso.y2.bin));
  lasso.pop.tot <- c(intercept=sum(acs13b$pwgtp), mu=sum(lasso.fit.pop*acs13b$pwgtp));

  if (sum(coef.lasso.y2.bin != 0)>1) {
    nhat = sum(wts); muhat = sum(wts*lasso.fit.samp); muhatsq = sum(wts*(lasso.fit.samp^2));
    Nhat = lasso.pop.tot[1]; MUhat = lasso.pop.tot[2];
    lambda <- as.vector(t(c(Nhat - nhat, MUhat - muhat))%*%solve(matrix(c(nhat, muhat, muhat, muhatsq),nrow=2)));
    g <- as.vector(1 + lambda[1] + lambda[2]*lasso.fit.samp);
    wts.new <- as.vector(wts*g);
    That.y2.bin=sum(wts.new*as.numeric(y2.bin.samp==1));
  } else {
    # intercept only model
    wts.new <- wts*lasso.pop.tot[2]/sum(wts*lasso.fit.samp);
    That.y2.bin=sum(wts.new*as.numeric(y2.bin.samp==1));
  }
  }, error=function(err) {
    possibleError=1;
  }
  );
  if (possibleError) {
#    print(sprintf("%s bootstrap error %s\n",sim.name,date()));
    return(c(NA));
  }

#  print(sprintf("%s complete 1 bootstrap %s\n",sim.name,date()));
  return(c(That.y2.bin));
}

