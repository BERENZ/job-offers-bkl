################################################################ 
# File:      functions.R
# Created:   12/3/2009
# Revised:   01/11/2015
# Author:    Jack Chen
#
# Description: Functions used by Simluation Project
################################################################ 
library(glmnet);
library(ROCR);
library(sampling);
library(survey);
#install.packages("devtools")
#devtools::install_github("tdhock/WeightedROC")
library(WeightedROC);
library(foreach);
library(doMC);

#library(Rcpp);
#sourceCpp("./mvrnormArma.cpp");

#cpprbinom <- cxxfunction(signature(ns="integer", ks="integer", ps="double"), plugin = "Rcpp", body='
#   int n = Rcpp::as<int>(ns);
#   int k = Rcpp::as<int>(ks);
#   double p = Rcpp::as<double>(ps);
#   Rcpp::RNGScope tmp;
#   Rcpp::NumericVector draws = Rcpp::floor(Rcpp::runif(n*k)+p);
#   return Rcpp::NumericMatrix(n, k, draws.begin());
#');
#
#cpprbinom <- cppFunction('NumericVector cpprbinom(int n, NumericVector p) {
#  NumericVector out(n);
#  NumericVector u = runif(n,0,1);
#  for (int i = 0; i < n; i++) {
#    out[i]=floor(u[i]+p[i]);
#  }
#  return out;
#}');

#expit2 <- cppFunction('NumericVector expit2(NumericVector x) {
#  NumericVector out(x.size());
#  for (int i = 0; i < x.size(); i++) {
#    out[i] = 1/(1+exp(-x[i]));
#  }
#  return out;
#}');

#library(rbenchmark);
#lala<-runif(1000);
#res <- benchmark(rbinom(1000,1,lala),cpprbinom(1000,lala));

# function      : ntile
# description   : break a numeric vector into bins
#
# args
# x             : a numeric vector
# n             : number of bins
#
# return value  : a vector of group numbers
ntile <-
function(x,n=4) {
  floor((n * (rank(x,ties.method = "first") - 1)/length(x)) + 1)
}

#  function     : propensity.select.opt (optimized version)
#  description  : select, based on "predilection" in reference
#
#  args
#  pop          : sample frame
#  pk           : unit response propensity in sample frame  
#  glm.obj      : glm object that modeled propensity to be in ref survey
#  n            : desired sample size
#  G            : number of classes
#
#  return value : a list containing vectors of selected sample index, one per class
#                 a vector containing number of selected sample in each class
#                 number of classes
#   
propensity.select2 <-
function(pop,pk,glm.obj,n,G=5) {
  p <- predict(glm.obj,type="response");
  breakpoints <- quantile(p, probs=seq(0,1,1/G));
  breakpoints[1]=0;
  breakpoints[G+1]=1;

  # propensity class sizes
  class.size <- rep(0,G);
  for (i in 1:G) {
    if (i != G) {
      class.size[i]=round(n/G);
    }
    else {
      # in case n is not divisible by G, last group has different size
      class.size[i]=n-sum(class.size[1:(G-1)]);
    }
  }

  # vector of selected sizes for each propensity class
  class.selected.size <- rep(0,G);

  # list of selected index for each propensity class
  class.index <- list(); class.index[[G+1]]=NA;

  # sample a bunch
  pop.select.ind <- which(UPpoisson(pk)==1);
  pop.select.pk <- predict(glm.obj,newdata=pop[pop.select.ind,],type="response");

  propensity.select.n = 0;
  for (i in 1:length(pop.select.ind)) {
    try.ind = pop.select.ind[i];
    try.pk = pop.select.pk[i];
    try.class = min(which(try.pk < breakpoints))-1;
#print(i);
#print(breakpoints);
#print(try.pk);
#print(try.class);
    if (class.selected.size[try.class] < class.size[try.class]) {
      class.index[[try.class]] <- c(class.index[[try.class]],try.ind);
      propensity.select.n=propensity.select.n+1;
    }
    class.selected.size[try.class]=class.selected.size[try.class]+1;
    if (propensity.select.n == n) break;
  }

  if (propensity.select.n != n) {
    stop("propensity.select: Not enough sample selected");
  }

  return(list(class.index=class.index, class.selected.size=class.selected.size));
}

#  function     : propensity.select
#  description  : select, based on "predilection" in reference
#
#  args
#  pop          : sample frame
#  pk           : unit reponse propensity in sample frame  
#  glm.obj      : glm object that modeled propensity to be in ref survey
#  n            : desired sample size
#  xtile        : vector of percentile breakpoints, default to quintile, (0,0.2,0.4,0.6,0.8,1)
#
#  return value : a list containing vectors of selected sample index, one per class
#                 a vector containing number of selected sample in each class
#                 number of classes
#   
propensity.select <-
function(pop,pk,glm.obj,n,xtile=NULL) {
  if (length(xtile)==0) {
    xtile=seq(0,1,0.2);
  }
  p <- predict(glm.obj,type="response");
  breakpoints <- quantile(p,xtile);

  # number of propensity classes
  G = length(xtile)-1;

  # propensity class sizes
  class.size <- rep(0,G);
  for (i in 1:G) {
    if (i != G) {
      class.size[i]=round(n/G);
    }
    else {
      # in case n is not divisible by G, last group has different size
      class.size[i]=n-sum(class.size[1:(G-1)]);
    }
  }

  # vector of selected sizes for each propensity class
  class.selected.size <- rep(0,G);

  # list of selected index for each propensity class
  class.index <- list(); class.index[[G+1]]=NA;

  propensity.select.n = 0;
  while(propensity.select.n < n) {
    try.ind = poi.fcn(pk,1);
    try.pk = predict(glm.obj,newdata=pop[try.ind,],type="response");
    try.class = min(which(try.pk <= breakpoints)) - 1;
    if (class.selected.size[try.class] < class.size[try.class]) {
      class.index[[try.class]] <- c(class.index[[try.class]],try.ind);
      propensity.select.n=propensity.select.n+1;
    }
    class.selected.size[try.class]=class.selected.size[try.class]+1;
  }

  return(list(class.index=class.index, class.selected.size=class.selected.size, G=G));
}

#  function     : formula.factor 
#  description  : construct formula of factored effects
#
#  args
#  v            : a vector of main effects
#
#  return value : a formula with each element in v assigned as factor
#   
formula.factor <-
function(v) {
  model <- NULL;
  for (i in 1:length(v)) {
    if (i == 1) {
      model <- paste0(model,"~factor(",v[i],")");
    }
    else {
      model <- paste0(model,"+factor(",v[i],")");
    }
  }
  return(formula(model));
}

#  function     : formula.str
#  description  : construct formula string
#
#  args
#  v            : a vector of variable names
#
#  return value : a string with for formula
#   
formula.str <-
function(v) {
  model.str <- NULL;
  for (i in 1:length(v)) {
    if (v[i] != "(Intercept)") {
      if (is.null(model.str)) {
        model.str <- paste0("~",v[i]);
      }
      else {
        model.str <- paste0(model.str,"+",v[i]);
      }
    }
  }
  return(model.str);
}

#  function     : expit
#  description  : inverse of logit
#
#  args
#  x            : a real number
#
#  return value : 1/(1+exp(-x))
#   
expit <-
function(x) {
  return(1/(1+exp(-x)));
}

#  function     : my.log
#  description  : log function to handle 0 value
#
#  args
#  x            : a real number
#
#  return value : log(.Machine$double.eps) if x is 0, log(x) if not
#   
my.log <-
function(x) {
  x[which(x == 0)] = .Machine$double.eps;
  return(log(x));
}

# function     :  draw.betas
# description  :  draw beta parameters from normal distributions
#
# args
# p            : number of betas
# beta.zero    : indices of betas that are preset to 0
# means        : vector of means; if single number,
#                then all betas share same mean
# sds          : vector of standard errors; if single number,
#                then all betas share same mean
# seed         : random seed
#
# return value : vector of betas drawn from independent normal distributions
#                from pre-specified means and stand deviations
#
draw.betas <-
function(p, beta.zero=NULL, means=NULL, sds=NULL, seed=0) {

  if (length(seed) == 1) {
    seed <- rep(seed, p);
  }

  # make sure not all betas are set to 0
  stopifnot(p > length(beta.zero));
  stopifnot(p == length(seed));

  # prespecified means for beta, default is 0
  if (is.null(means)) {
    means <- rep(0,p);
  } else {
    if (length(means) > 1) {
      stopifnot(length(means) == p);
    } else {
      means <- rep(means,p);
    }
  }

  # prespecified sd for beta, default is 1
  if (is.null(sds)) {
    sds <- rep(1,p);
  } else {
    if (length(sds) > 1) {
      stopifnot(length(sds) == p);
    } else {
      sds <- rep(sds,p);
    }
  }

  # draw independent normals
  betas <- rep(0, p);
  for (i in 1:p) {
    if (!(i %in% beta.zero)) {
      set.seed(seed[i]+i);
      betas[i] <- rnorm(1, mean=means[i], sd=sds[i]);
    }
  }

  return(betas);
}

# function     : draw.ys
# description  : draw independent bernoulli random variables
#
# args
# ps           : vector of bernoulli probabilities
# seed         : random seed
#
# return value : vector of independent bernoulli random variables
#
draw.ys <-
function(ps, seed=0) {
  set.seed(seed);
  ys <- rep(NA, length(ps));
  for (i in 1:length(ps)) {
    ys[i] = rbinom(1,1,ps[i]);
  }
  return(ys);
}

# function     : generate.x
# description  : generate X with variables
#                X1 ~ unif(X1.min, X1.max);
#                X2 = exp(
#                X3 ~ N(X3.mean, X3.sd^2)
#                X4 ~ X3^2
#                C ~ multinom(k/sum(k)); k = 1, ..., C.num.cat
#
# args
# n            : number of data points to generate
# X1.min       : min of uniform X1
# X1.max       : max of uniform X1
# X2.lambda    : lambda of exponential X2
# X3.mean      : mean of normal X3, default is 0
# X3.sd        : standard deviation of normal X3, default is 1
# C.cat        : vector of number of categories in categorical variables
#                each category has probability of 1 = category #/sum(category #)
# seed         : random seed
#
# return value : data matrix with X1, X2, X3, C1, and C2 as described
#
generate.x <-
function(n, X1.min=1, X1.max=3, X2.lambda=1, X3.mean=0, X3.sd=1, C.cat=2, seed=0) {
  set.seed(seed);

  # continuous uniform
  X1 <- runif(n, X1.min, X1.max);

  # continuous exponential
  X2 <- rexp(n, X2.lambda);

  # continuous normal
  X3 <- rnorm(n, mean=X3.mean, sd=X3.sd);

  # categorical
  C.cat <- as.vector(C.cat);
  C.length <- length(C.cat);
  C.matrix <- as.data.frame(matrix(rep(NA, n*C.length), nrow=n));
  for (i in 1:length(C.cat)) {
    set.seed(seed + i);
    Ci.dummy <- t(rmultinom(n, 1, c(1:C.cat[i])/sum(1:C.cat[i])));
    Ci <- as.factor(Ci.dummy%*%c(1:C.cat[i]));
    C.matrix[,i] <- Ci;
  }
  dimnames(C.matrix)[[2]] <- paste("C", c(1:C.length), sep="");

  X <- data.frame(X1,X2,X3,C.matrix);
  return(X);
}

# function     : generate.noise.unif
# description  : generate p uniform random vectors
#
# args
# n            : number of data points
# p            : number of uniform distribution
# u.min,u.max  : lower/upper bounds of uniform distribution
#                can be single number for all p distributions
#                or vector of numbers for each of the distribution
# colname      : prefix of column names in final data matrix
# seed         : random seed
#
# return value : data matrix of dimension (n x p). Each column corresponds to
#                a uniform distribution with specified min/max
#
generate.noise.unif <-
function(n, p=2, u.min=0, u.max=1, seed=c(1,2), colname="U.cont"){

  stopifnot(length(seed) == p);

  # make sure min and max of the distribution(s)
  # is specified correctly
  stopifnot(length(u.min) == length(u.max));

  if (length(u.min) == 1) {
    u.min <- rep(u.min, p);
    u.max <- rep(u.max, p);
  }
  stopifnot(sum(u.min < u.max) == p);

  U <- as.data.frame(matrix(rep(NA, p*n), nrow=n));
  for (i in 1:p) {
    set.seed(i);
    U[,i] <- runif(n, u.min[i], u.max[i]);
  }

  dimnames(U)[[2]] <- paste(colname, c(1:p), sep="");
  return(U);
}

# function     : generate.noise.multi
# description  : generate p multinomial random variables
#
# args
# n            : number of data points
# p            : number of categorical variables, each with
#                a multinomial distribution of 1/num.cat
# num.cat      : number of categories
#                can be single number for all p distributions
#                or vector of numbers for each of the distribution
# colname      : prefix of column names in final data matrix
# seed         : random seed
#
# return value : data matrix of dimension (n x p). Each column corresponds to
#                a multinomial distribution with probabilities 1/num.cat
#
generate.noise.multi <-
function(n, p=2, num.cat=5, seed=0, colname="U.cat"){

  if (length(seed) == 1) {
    seed <- rep(seed, p);
  }
  stopifnot(length(seed) == p);

  if (length(num.cat) == 1) {
    num.cat <- rep(num.cat, p);
  }
  stopifnot(length(num.cat) == p);

  U <- as.data.frame(matrix(rep(NA, p*n), nrow=n));
  for (i in 1:p) {
    set.seed(seed[i]);
    Ui.dummy <- t(rmultinom(n, 1, rep(1/num.cat[i], num.cat[i])));
    U[,i] <- as.factor(Ui.dummy%*%c(1:num.cat[i]));
  }

  dimnames(U)[[2]] <- paste(colname, c(1:p), sep="");
  return(U);
}


# function     : generate.cor.unif
# description  : generate p uniform distributions correlated to X
#
# args
# X            : a vector of normal values, or matrix of normal vectors
# rho          : desired Pearson correlation, must be length p.
# seed         : random seed
#
# return value : data matrix of dimension same as X. Each column corresponds to
#                a uniform vector correlated to corresponding column in X
#                with correlation specified by rho
#
# ref: http://comisef.wikidot.com/tutorial:correlateduniformvariates
#
generate.cor.unif <-
function(X, rho, seed=0, colname="P.unif"){

  X <- as.matrix(X);
  n = nrow(X);
  p = ncol(X);

  # input error check
  if (length(seed) == 1) {
    seed <- rep(seed, p);
  }

  stopifnot(length(rho) == p);
  stopifnot(length(seed) == p);

  P <- as.data.frame(matrix(rep(NA, n*p), nrow=n));

  rho.adj <- c(1, 2*sin(pi*rho/6));

  # create correlation matrix

  # the correlation of original normals
  M.X <- cor(X, X);

  # the overall correlation matrix where first p vectors
  # or the original
  M <- matrix(rep(0, 4*p^2), nrow=2*p);
  M[1:p,1:p] <- M.X;

  # the (p+i)th column in C has correlation rho[i] with (i)th
  # column in C; no correlations between uniforms
  for (i in 1:p) {
    M[i, p+i] = 2*sin(pi*rho[i]/6);
  }

  # transformation adjustment (see reference)
  M <- 2*sin(pi * M / 6)
  diag(M) <- rep(1, 2*p);

  # normalize random normals in X,
  # generate p random normals to go along with the p specified in X
  Z <- matrix(rep(NA, n*2*p), nrow=n);

  for (i in 1:p) {
    mean.Xi <- mean(X[,i]);
    sd.Xi <- sd(X[,i]);
    Z[,i] <- (X[,i] - mean.Xi)/sd.Xi;
  }

  # generate random normals for new vectors
  for (i in 1:p) {
    set.seed(seed[i]);
    Z[,p+i] <- rnorm(n)
  }

  C <- chol(M);
  Y <- Z %*% C;

  P <- as.data.frame(pnorm(Y[, (p+1):(2*p)]));

  dimnames(P)[[2]] <- paste(colname, c(1:p), sep="");
  return(P);
}

# function     : generate.cor.norm
# description  : generate p normal distributions correlated to X
#
# args
# X            : a vector of normal values, or matrix of
#                p normal vectors
# rho          : correlation with X. Length must equal to
#                the number of normal vectors in X
# seed         : random seed
#
# return value : data matrix of dimension same as X. Each column corresponds to
#                a normal vector correlated to corresponding column in X
#                with correlation specified by rho
#
# ref:
#      Goldman, R.N., and McKenzie, J.D. (2009),
#      "Creating Realistic Data Sets with Specified Properties via Simulation",
#      Teaching Statistics, 31, 1, 7-10
#
generate.cor.norm <-
function(X, rho, seed=0, colname="P.norm"){

  X <- as.matrix(X);
  n = nrow(X);
  p = ncol(X);

  # seed used to generate each normal
  if (length(seed) == 1) {
    seed <- rep(seed, p);
  }

  # input error check
  stopifnot(length(rho) == p);
  stopifnot(length(seed) == p);

  # result matrix
  P <- as.data.frame(matrix(rep(NA, n*p), nrow=n));

  for (i in 1:p) {
    mean.Xi = mean(X[,i]);
    sd.Xi = sd(X[,i]);

    # normalize X
    Z1 <- (X[,i] - mean.Xi)/sd.Xi;

    set.seed(seed[i]);
    Z2 <- rnorm(n);

    # create correlated normals
    # Z1 and Z3 have rho[i] correlation
    Z3 = rho[i]*Z1 + sqrt((1-rho[i]^2))*Z2;

    # transform back to get same scale
    P[,i] <- Z3*sd.Xi + mean.Xi;
  }

  dimnames(P)[[2]] <- paste(colname, c(1:p), sep="");
  return(P);
}

# function     : generate.permute
# description  : permute the entries of each column in X
#
# args
# X            : a vector, or matrix
# prop         : proportion of elements of each column in X to be permuted
# seed         : random seed
#
# return value : data matrix of dimension same as X. Each column corresponds to
#                a factor with a portion permuted
#
generate.permute <-
function(X, prop, seed=0, colname="P.perm"){

  X <- as.matrix(X);
  n = nrow(X);
  p = ncol(X);

  # input error check
  if (length(seed) == 1) {
    seed <- rep(seed, p);
  }

  stopifnot(length(prop) == p);
  stopifnot(length(seed) == p);

  P <- as.data.frame(matrix(rep(NA, n*p), nrow=n));

  for (i in 1:p) {
    set.seed(seed[i]+i*seed[i]);
    change.ind <- sample(1:n, n*prop[i]);

    set.seed(seed[i]+i);
    change.to.ind <- sample(1:n, n*prop[i]);
    P[,i] <- X[,i];
    P[change.ind,i] <- X[change.to.ind,i];
  }

  dimnames(P)[[2]] <- paste(colname, c(1:p), sep="");
  return(P);
}

# function     : rmse
# description  : calculate empirical bias, var, and rmse
#                from simulation outputs
#
# args
# results      : simulation outputs
# T.pop        : True population quantity
#
# return value : a list of bias, variance, and rmse
#
rmse <-
function(results, T.pop) {
  num.na = sum(is.na(results));
  v <- var(results, na.rm=T);
  bias <- mean(results, na.rm=T) - T.pop;
  return(list(bias=bias, v=v, rmse=sqrt(bias^2+v), num.na=num.na));
}

# function     : pps.fcn
# description  : Select pps sample of size n using Hartley-Rao algorithm.
#                Returns indices in population of sample units.
#                If any units are certainties, then NULL vector is returned.
#
# args
# x            : population
# n            : sample size
#
# return value : a vector of indices
#
pps.fcn <-
function(x, n) {
    N <- length(x);
    if (n > N) {
        stop("Sample size > pop.size.");
    }
    cumsums <- cumsum(x);
    Skip <- cumsums[N]/n;
    if (max(x) > Skip) {
        cat("At least 1 unit is a certainty. Sample not selected.\n");
        indices <- NULL;
    } else {
        R <- runif(1, 0, Skip);
        u <- R + Skip*(0:(n-1));
        indices <- N + 1 - outer(u, cumsums, "<=") %*%
          matrix(rep(1,N), ncol = 1);
    }
    return(indices);
}

# function     : pps.random.fcn
# description  : Select pps sample of size n using Hartley-Rao algorithm after
#                randomizing population order.  Returns indices in population
#                of sample units.
#
# args
# x            : population vector of sizes
# n            : sample size
#
# return value : vector of selected indices
#
pps.random.fcn <-
function(x, n) {
    N <- length(x);
    if (n > N) {
        stop("Sample size > pop size.");
    }
    X <- cbind(x, 1:N);
    X.tmp <- X[sort.list(sample(1:N, N)),];
    smp.tmp <- pps.fcn(X.tmp[,1], n);
    ind <- X.tmp[smp.tmp, 2];
    return(sort(ind));
}

# function     : poi.fcn
# description  : Select poisson sample of size n
#
# args
# pk           : vector of probabilities having desired characteristics
#                e.g. probabilities of being a volunteer web respondent
# n            : desired sample size (fixed, by inducing some dependency)
#
# return value : vector of selected indices
#
poi.fcn <-
function(pk, n) {
    En = sum(pk);  # expected number of cases with desired characteristics
    if (n > length(pk)) {
        stop("Sample size > pop size.");
    }

    n.remaining = n;
    ind <- c();
    
    while (n.remaining > 0) {
      try.ind <- which(UPpoisson(n.remaining*pk/En)==1);
      ind <- union(ind,try.ind);
      n.remaining = n - length(ind);
    }

    # in case more than n selected, choose first n
    return(sort(ind[1:n]));
}

# function     : poi.fcn.rnd
# description  : Select poisson sample of size n
#
# args
# pk           : vector of probabilities having desired characteristics
#                e.g. probabilities of being a volunteer web respondent
# n            : desired sample size (random)
#
# return value : vector of selected indices
#
poi.fcn.rnd <-
function(pk, n) {
    En = sum(pk);  # expected number of cases with desired characteristics
    if (n > length(pk)) {
        stop("Sample size > pop size.");
    }

    ind <- which(UPpoisson(n.remaining*pk/En)==1);
    return(sort(ind));
}

# returns a data frame of two variables which correlate with a population correlation of rho
# If desired, one of both variables can be fixed to an existing variable by specifying x
getBiCop <- function(n, rho, mar.fun=rnorm, x = NULL, ...) {
     if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(n, ...)}
     if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")

     C <- matrix(rho, nrow = 2, ncol = 2)
     diag(C) <- 1

     C <- chol(C)

     X2 <- mar.fun(n)
     X <- cbind(X1,X2)

     # induce correlation (does not change X1)
     df <- X %*% C

     ## if desired: check results
     #all.equal(X1,X[,1])
     #cor(X)

     return(df)
}

# n      sample size
# rho    desired correlation = cos(angle)
# x      fixed given data
generate.cor <- function(n, rho, x, mar.fun=rnorm) {
  if (!is.null(x)) {x1 <- x} else {x1 <- mar.fun(n, ...)}
  if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
  theta <- acos(rho)             # corresponding angle
  x2    <- mar.fun(n);
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

#  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
#  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
#  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
#  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
#  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

  x2o <- rep(NA,n);
  for (i in 1:n) {
    print(i);
    QQt_i <- Q[i] * Q # t cross product element
    QQt_i[i] = 1 - QQt_i[i];
    x2o[i] = sum(QQt_i * Xctr[, 2]);
  }

  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  cor(x1, x)                                    # check correlation = rho
  return(x);
}

# function       : wabserr
# description    : weighted absolute error

# function       : mywauc
# description    : weighted Area Under Curve
mywauc <- function (prob, y, w) 
{
    if (missing(w)) {
        rprob = rank(prob)
        n1 = sum(y)
        n0 = length(y) - n1
        u = sum(rprob[y == 1]) - n1 * (n1 + 1)/2
        u/(n1 * n0)
    }
    else {
        rprob = runif(length(prob))
        op = order(prob, rprob)
        y = y[op]
        w = w[op]
        cw = cumsum(w)
        w1 = w[y == 1]
        cw1 = cumsum(w1)
        wauc = sum(w1 * (cw[y == 1] - cw1))
        sumw = cw1[length(cw1)]
        sumw = sumw * (cw[length(cw)] - sumw)
        wauc/sumw
    }
}

# function       : wabserr
# description    : weighted absolute error
wabserr <- function(yhat, y, w) {
  sum(abs(yhat - y)*w)/sum(w);
};

# function       : mycv.glmnet
# description    : performs cross-validation for lasso regression
#
# args
# x              : candidate variables
# y              : dependent variable
# standardize    : whether to standardize (subtract mean and divide by sd, FALSE if all variables are categorical)
# intercept      : whether to include intercept in model
# penalty.factor : penalty factor associated with each candidate variable (default to 1)
# seed           : seed used to generate cross-validation samples
# nfolds         : number of cross-validation groups
# alpha          : penalty coefficient, 0=bridge regression, 1=lasso
# family         : distribution family (binomial or gaussian supported)
# fun            : function to determine the best metric (max for AUC, min for rmse)
# type           : predicted value type (response = predicted mean)
# nlambda        : number of lambda grids
# weights        : observation weights
# cores          : number of cpu cores for parallel processing
# adaptive       : use adaptive Lasso (default=FALSE)
#
# return value : coef =lasso coefficient corresponding to optimal penalty parameter obtained through cross-validation
#                metrics = matrix of metrics, each row is a lambda, each column is measure on a test set (# of column = # of CV)
#                lambda = vector of lambda grids
#                best.metric.index = the index (row number in metrics) the has the best average measure
#
# assume x does not contain intercept column, and all variables are dummy vectors (glmnet does not dummify variables internally)
mycv.glmnet <-
function (x, y, standardize = FALSE, intercept = TRUE, penalty.factor = NULL, 
    seed = NULL, nfolds = 5, family = "binomial", alpha = 1, 
    measure = "auc", fun = "max", type = "response", weights = NULL, 
    nlambda = 100, cores = 1, adaptive = FALSE) {

    if (family == "binomial") {
      y <- as.factor(y);
    }
    if (length(penalty.factor) == 0) {
      penalty.factor = rep(1, ncol(x));
    }
    if (is.null(weights)) {
      weights = rep(1, nrow(x));
    }
    if (!is.null(seed)) set.seed(seed);

    N = length(y);
    foldsize = round(N/nfolds);
    fold.index <- list();
    fold.index[[nfolds]] <- NULL;
    index <- c(1:N);
    for (i in 1:(nfolds - 1)) {
      if (i == 1) {
        samp.index <- sample(index, size = foldsize, replace = FALSE);
      } else {
        samp.index <- sample(index[-unlist(fold.index)], size = foldsize, replace = FALSE);
      }
      fold.index[[i]] <- samp.index;
    }
    fold.index[[nfolds]] = index[-unlist(fold.index)];

    glmnet.obj0 <- NULL;
    if (adaptive == FALSE) {
      glmnet.obj0 <- glmnet(x, y, penalty.factor = penalty.factor, 
          family = family, alpha = alpha, standardize = standardize, 
          intercept = intercept, weights = weights, nlambda = nlambda);

      metrics0 <- matrix(rep(NA, nfolds * length(glmnet.obj0$lambda)), ncol = nfolds);

      for (i in 1:nfolds) {
        train.index <- unlist(fold.index[-i]);
        test.index <- fold.index[[i]];
        train.glmnet <- glmnet(x[train.index, ], y[train.index], 
                penalty.factor = penalty.factor, family = family, 
                lambda = glmnet.obj0$lambda, standardize = standardize, 
                intercept = intercept, alpha = alpha,
                weights = weights[train.index], nlambda = nlambda);
        predict.y <- predict(train.glmnet, x[test.index, ], type = type);
        metric.fold <- rep(NA, length(glmnet.obj0$lambda));

        for (j in 1:length(glmnet.obj0$lambda)) {
          if (family == "binomial") {
#            tp.fp <- WeightedROC(predict.y[, j], y[test.index], weights[test.index]);
#            wauc <- WeightedAUC(tp.fp);
            metric.fold[j] <- glmnet::auc(y[test.index],predict.y[,j],weights[test.index]);
#            metric.fold[j] <- eval(call(fun, wauc));
          } else {
            metric.fold[j] <- sum(abs(predict.y[, j] - y[test.index]) *
                (weights[test.index]))/sum(weights[test.index]);
          }
        }
        metrics0[, i] <- metric.fold
      }

      metric0.mean <- apply(metrics0, 1, mean);
      best.metric0 = eval(call(fun, metric0.mean));
      best.metric0.index = which(metric0.mean == best.metric0)[1];
      return(list(coef = coef(glmnet.obj0)[, best.metric0.index], 
          metrics = metrics0, lambda = glmnet.obj0$lambda, 
          best.metric.index = best.metric0.index, method = 0))
    }

    predictors <- colnames(x);
    glm.formula <- update(formula(formula.str(predictors)),"y~.");

    # adaptive
    glmnet.obj1 <- NULL
    glmnet.obj2 <- NULL
    glmnet.obj3 <- NULL
    glmnet.obj4 <- NULL

    taus <- c(0.1,0.5,1,2);
    tau1 = taus[1];
    tau2 = taus[2];
    tau3 = taus[3];
    tau4 = taus[4];

    # initial weights, whole data
    data.ls.fit <- data.frame(y=y, x);
    ls.fit <- NULL;
    if (family == "binomial") {
      dsgn <- svydesign(ids=~1,weights=~weights,data=data.ls.fit);
      ls.fit <- svyglm(glm.formula, design=dsgn, family=quasibinomial());
    } else {
      ls.fit <- lm(glm.formula, data = data.ls.fit, weights = weights);
    }
    beta.init.whole <- (coef(ls.fit))[-1];

    # initial weights, each fold
    beta.list <- list(); beta.list[[nfolds]] <- NULL;
    for (i in 1:nfolds) {
      train.index <- unlist(fold.index[-i]);
      test.index <- fold.index[[i]];
      data.ls.fit <- data.frame(y=y[train.index], x[train.index,]);
      ls.fit <- NULL
      if (family == "binomial") {
        dsgn <- svydesign(ids=~1,weights=~weights[train.index],data=data.ls.fit);
        ls.fit <- svyglm(glm.formula, design=dsgn, family=quasibinomial());
      } else {
        ls.fit <- lm(glm.formula, data = data.ls.fit, weights = weights[train.index]);
      }
      beta.list[[i]] <- coef(ls.fit)[-1];
    }

    glmnet.objs <- list();
#    foreach (i=1:4) %dopar% {
    for (i in 1:4) {
      glmnet.objs[[i]] <- glmnet(x, y, penalty.factor = (1/abs(beta.init.whole))^taus[i], 
              family = family, alpha = alpha, standardize = standardize, 
              intercept = intercept, weights = weights, nlambda = nlambda);
    }
    glmnet.obj1 <- glmnet.objs[[1]];
    glmnet.obj2 <- glmnet.objs[[2]];
    glmnet.obj3 <- glmnet.objs[[3]];
    glmnet.obj4 <- glmnet.objs[[4]];


#    metrics1 <- matrix(rep(NA, nfolds * length(glmnet.obj1$lambda)), ncol = nfolds);
#    metrics2 <- matrix(rep(NA, nfolds * length(glmnet.obj2$lambda)), ncol = nfolds);
#    metrics3 <- matrix(rep(NA, nfolds * length(glmnet.obj3$lambda)), ncol = nfolds);
#    metrics4 <- matrix(rep(NA, nfolds * length(glmnet.obj4$lambda)), ncol = nfolds);

    metrics <- list();
    metrics[[4]] <- NULL;
    for (k in 1:4) {

      metric.folds <- matrix(rep(NA, nfolds*length(glmnet.objs[[k]]$lambda)), ncol=nfolds);
      for (i in 1:nfolds) {
        train.index <- unlist(fold.index[-i]);
        test.index <- fold.index[[i]];
        beta.init <- beta.list[[i]];
        train.glmnet <- glmnet(x[train.index, ], y[train.index], 
                penalty.factor = (1/abs(beta.init))^taus[k], family = family, 
                lambda = glmnet.objs[[k]]$lambda, standardize = standardize, 
                intercept = intercept, alpha = alpha,
                weights = weights[train.index], nlambda = nlambda);
        predict.y <- predict(train.glmnet, x[test.index, ], type = type);
        metric.fold <- rep(NA, length(glmnet.objs[[k]]$lambda));
        grid.n <- min(length(train.glmnet$lambda), length(glmnet.objs[[k]]$lambda));
        for (j in 1:grid.n) {
          if (family == "binomial") {
            metric.fold[j] <- glmnet::auc(y[test.index],predict.y[,j],weights[test.index]);
          } else {
            metric.fold[j] <- sum(abs(predict.y[, j] - y[test.index]) *
                (weights[test.index]))/sum(weights[test.index]);
          }
        }
        metric.folds[, i] <- metric.fold
      }
      metrics[[k]] <- metric.folds;
    }
    metrics1 <- metrics[[1]];
    metrics2 <- metrics[[2]];
    metrics3 <- metrics[[3]];
    metrics4 <- metrics[[4]];

    metric1.mean <- apply(metrics1, 1, mean, na.rm = T);
    metric2.mean <- apply(metrics2, 1, mean, na.rm = T);
    metric3.mean <- apply(metrics3, 1, mean, na.rm = T);
    metric4.mean <- apply(metrics4, 1, mean, na.rm = T);
    best.metric1 = eval(call(fun, metric1.mean));
    best.metric2 = eval(call(fun, metric2.mean));
    best.metric3 = eval(call(fun, metric3.mean));
    best.metric4 = eval(call(fun, metric4.mean));
    best.metric = eval(call(fun, best.metric1, best.metric2, best.metric3, best.metric4));
    best.metric.method = which(c(best.metric1, best.metric2, best.metric3, best.metric4) == best.metric)[1];
    if (best.metric.method == 1) {
      best.metric1.index = which(metric1.mean == best.metric1)[1];
      return(list(coef = coef(glmnet.obj1)[, best.metric1.index], 
                metrics = metrics1, lambda = glmnet.obj1$lambda, 
                best.metric.index = best.metric1.index, method = 1));
    }
    if (best.metric.method == 2) {
      best.metric2.index = which(metric2.mean == best.metric2)[1]
      return(list(coef = coef(glmnet.obj2)[, best.metric2.index], 
                metrics = metrics2, lambda = glmnet.obj2$lambda, 
                best.metric.index = best.metric2.index, method = 2));
    }
    if (best.metric.method == 3) {
      best.metric3.index = which(metric3.mean == best.metric3)[1]
      return(list(coef = coef(glmnet.obj3)[, best.metric3.index], 
                metrics = metrics3, lambda = glmnet.obj3$lambda, 
                best.metric.index = best.metric3.index, method = 3));
    }
    best.metric4.index = which(metric4.mean == best.metric4)[1]
    return(list(coef = coef(glmnet.obj4)[, best.metric4.index], 
                metrics = metrics4, lambda = glmnet.obj4$lambda, 
                best.metric.index = best.metric4.index, method = 4));
}

# specify method in adaptive
mycv.glmnet.method <-
function (x, y, standardize = FALSE, intercept = TRUE, penalty.factor = NULL, 
    seed = NULL, nfolds = 5, family = "binomial", alpha = 1, 
    measure = "auc", fun = "max", type = "response", weights = NULL, 
    nlambda = 100, cores = 1, method=1) {

#    x.mat <- as.matrix(x);
    predictors <- colnames(x);
    glm.formula <- update(formula(formula.str(predictors)),"y~.");

    if (family == "binomial") {
        y <- as.factor(y);
    }
    if (length(penalty.factor) == 0) {
        penalty.factor = rep(1, ncol(x));
    }
    if (is.null(weights)) {
        weights = rep(1, nrow(x));
    }
    if (!is.null(seed)) set.seed(seed);

    N = length(y);

    foldsize = round(N/nfolds);
    fold.index <- list();
    fold.index[[nfolds]] <- NULL;
    index <- c(1:N);
    for (i in 1:(nfolds - 1)) {
      if (i == 1) {
        samp.index <- sample(index, size = foldsize, replace = FALSE);
      } else {
        samp.index <- sample(index[-unlist(fold.index)], size = foldsize, replace = FALSE);
      }
      fold.index[[i]] <- samp.index;
    }
    fold.index[[nfolds]] = index[-unlist(fold.index)];

    # initial weights, whole data
    data.ls.fit <- data.frame(y=y, x);
    ls.fit <- NULL;
    if (family == "binomial") {
      dsgn <- svydesign(ids=~1,weights=~weights,data=data.ls.fit);
      ls.fit <- svyglm(glm.formula, design=dsgn, family=quasibinomial());
    } else {
      ls.fit <- lm(glm.formula, data = data.ls.fit, weights = weights);
    }
    beta.init.whole <- (coef(ls.fit))[-1];

    # initial weights, each fold
    beta.list <- list(); beta.list[[nfolds]] <- NULL;
    for (i in 1:nfolds) {
      train.index <- unlist(fold.index[-i]);
      test.index <- fold.index[[i]];
      data.ls.fit <- data.frame(y=y[train.index], x[train.index,]);
      ls.fit <- NULL
      if (family == "binomial") {
        dsgn <- svydesign(ids=~1,weights=~weights[train.index],data=data.ls.fit);
        ls.fit <- svyglm(glm.formula, design=dsgn, family=quasibinomial());
      } else {
        ls.fit <- lm(glm.formula, data = data.ls.fit, weights = weights[train.index]);
      }
      beta.list[[i]] <- coef(ls.fit)[-1];
    }

    taus <- c(0.1,0.5,1,2);

    glmnet.obj <- glmnet(x, y, penalty.factor = (1/abs(beta.init.whole))^taus[method], 
            family = family, alpha = alpha, standardize = standardize, 
            intercept = intercept, weights = weights, nlambda = nlambda);

    metrics <- matrix(rep(NA, nfolds*length(glmnet.obj$lambda)), ncol=nfolds);
    for (i in 1:nfolds) {
      train.index <- unlist(fold.index[-i]);
      test.index <- fold.index[[i]];
      beta.init <- beta.list[[i]];
      train.glmnet <- glmnet(x[train.index, ], y[train.index], 
              penalty.factor = (1/abs(beta.init))^taus[method], family = family, 
              lambda = glmnet.obj$lambda, standardize = standardize, 
              intercept = intercept, alpha = alpha,
              weights = weights[train.index], nlambda = nlambda);
      predict.y <- predict(train.glmnet, x[test.index, ], type = type);
      metric.fold <- rep(NA, length(glmnet.obj$lambda));
      grid.n <- min(length(train.glmnet$lambda), length(glmnet.obj$lambda));
      for (j in 1:grid.n) {
        if (family == "binomial") {
          metric.fold[j] <- glmnet::auc(y[test.index],predict.y[,j],weights[test.index]);
        } else {
          metric.fold[j] <- sum(abs(predict.y[, j] - y[test.index]) *
              (weights[test.index]))/sum(weights[test.index]);
        }
      }
      metrics[, i] <- metric.fold
    }

    metric.mean <- apply(metrics, 1, mean, na.rm = T);
    best.metric = eval(call(fun, metric.mean));
    best.metric.index = which(metric.mean == best.metric)[1];
    return(list(coef = coef(glmnet.obj)[, best.metric.index], 
                metrics = metrics, lambda = glmnet.obj$lambda, 
                best.metric.index = best.metric.index, method = method));
}

