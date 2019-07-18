library(foreign);
library(survey);

source("C:/Penndata/survey/Chen/NHIS/functions.R");


# Load NHIS13 data
nhis13 <- read.dta("C:/Penndata/survey/Chen/NHIS/process2.dta");
nhis13 <- nhis13[which(!is.na(nhis13$educ2)),];

nhis13[which(nhis13$age_p >= 18 & nhis13$age_p <= 29),"agegrp"]=1;
nhis13[which(nhis13$age_p >= 30 & nhis13$age_p <= 39),"agegrp"]=2;
nhis13[which(nhis13$age_p >= 40 & nhis13$age_p <= 49),"agegrp"]=3;
nhis13[which(nhis13$age_p >= 50 & nhis13$age_p <= 59),"agegrp"]=4;
nhis13[which(nhis13$age_p >= 60 & nhis13$age_p <= 69),"agegrp"]=5;
nhis13[which(nhis13$age_p >= 70 & nhis13$age_p <= 79),"agegrp"]=6;
nhis13[which(nhis13$age_p >= 80),"agegrp"]=7;
nhis13[which(nhis13$race == 4),"race"]=6;

nhis13$gender=nhis13$sex;
nhis13[which(is.na(nhis13$marst)),"marst"]=3;
nhis13$y2 <- nhis13$cancer;

nhis13$faminc_q <- nhis13$faminc_q2;
nhis13[which(is.na(nhis13$faminc_q2)),"faminc_q"] = 9;
#nhis13 <- subset(nhis13,faminc_q != 9);

nhis13 <- nhis13[!is.na(nhis13$y2),];
nhis13 <- nhis13[!is.na(nhis13$employed),];

#
#model <- formula(~factor(region)+factor(agegrp)+factor(gender)+factor(educ2)+factor(race)+factor(faminc_q));
#model <- formula(~factor(region)+factor(agegrp)+factor(gender)+factor(educ2)+factor(race)+factor(marst)+factor(employed)+factor(faminc_q));
#model <- formula(~factor(region)+factor(agegrp)+factor(gender)+factor(educ2)+factor(race)+factor(marst)+factor(employed)+factor(gender):factor(race)+factor(gender):factor(agegrp)+factor(agegrp):factor(race));
model <- formula(~factor(region)+factor(agegrp)+factor(gender)+factor(educ2)+factor(race)+factor(marst)+factor(employed)+factor(faminc_q))

model.GREG <- formula(~factor(agegrp)+factor(educ2)+factor(race)+factor(employed)+factor(faminc_q));

pop.bench.tot <- apply(model.matrix(model.GREG,acs13b)*acs13b$pwgtp,2,sum);
pop.bench.tot.full <- apply(model.matrix(model,acs13b)*acs13b$pwgtp,2,sum);

nhis13$wts2 <- nhis13$wtfa_sa * sum(acs13b$pwgtp)/sum(nhis13$wtfa_sa);

N = sum(acs13b$pwgtp);
n = nrow(nhis13);
wts = rep(N/n,n);
nhis13$wts <- wts;

# original design
dsgn2 <- svydesign(ids=~psu_p, strata=~strat_p, data=nhis13, weights=~wts2, nest=TRUE);
svytotal(~y2,dsgn2);

# HT assume SRS
dsgn2b <- svydesign(ids=~0, strata=NULL, data=nhis13, weights=~wts);
svytotal(~y2,dsgn2b);

# GREG weights

dsgn <- svydesign(ids=~0, strata=NULL, data=nhis13, weights=~wts);
sam.lin <- calibrate(design=dsgn, formula = model.GREG, population=pop.bench.tot,bounds=c(-Inf,Inf),calfun=c("linear"));
nhis13$wts3 <- weights(sam.lin);

dsgn3 <- svydesign(ids=~0, strata=NULL, data=nhis13, weights=~wts3);
svytotal(~y2,dsgn3);

sam.lin.full <- calibrate(design=dsgn, formula = model, population=pop.bench.tot.full,bounds=c(-Inf,Inf),calfun=c("linear"));
nhis13$wts4 <- weights(sam.lin.full);

dsgn4 <- svydesign(ids=~0, strata=NULL, data=nhis13, weights=~wts4);
svytotal(~y2,dsgn4);


# LASSO weights
X.y2.bench.lasso <- model.matrix(model,acs13b);

X.y2.lasso.samp <- model.matrix(model,nhis13);
wts.bench <- acs13b$pwgtp;
y2.samp <- nhis13$y2;

cvobj.y2.lasso <- mycv.glmnet(X.y2.lasso.samp[,-1],y2.samp,family="binomial",seed=2341,nfolds=5,alpha=1,measure="auc",type="response",fun="max",adaptive=FALSE,nlambda=100);

coef.y2.lasso <- cvobj.y2.lasso$coef;

lasso.fit.pop <- expit(as.vector(cbind(X.y2.bench.lasso) %*% as.vector(coef.y2.lasso)));
lasso.fit.samp <- expit(as.vector(cbind(X.y2.lasso.samp) %*% as.vector(coef.y2.lasso)));
lasso.pop.tot <- c(intercept=sum(wts.bench), mu=sum(lasso.fit.pop*wts.bench));

nhat = sum(wts); muhat = sum(wts*lasso.fit.samp); muhatsq = sum(wts*(lasso.fit.samp^2));
Nhat = lasso.pop.tot[1]; MUhat = lasso.pop.tot[2];
lambda <- as.vector(t(c(Nhat - nhat, MUhat - muhat))%*%solve(matrix(c(nhat, muhat, muhat, muhatsq),nrow=2)));
wts4 <- as.vector(wts*(1 + lambda[1] + lambda[2]*lasso.fit.samp));

sum(nhis13$y2* wts4);

data.samp <- nhis13;
data.bench <- acs13b;

data.samp$pk <- n/N;
pk.samp <- rep(n/N,n);

source("C:/Penndata/survey/Chen/NHIS/T.y2.lasso.R");

n.boot=500;
date();
  T.y2.lasso.boot <- boot(data=data.samp, statistic=T.y2.lasso, R=n.boot, ncpus=16, parallel="multicore",weights=1/pk.samp,sim.name="nhis13",alpha=1,nlambda=100,nfolds=5,adaptive=FALSE,method=0);
date();
v = var(T.y2.lasso.boot$t[,1],na.rm=TRUE);
m = mean(T.y2.lasso.boot$t[,1],na.rm=TRUE);
q <- quantile(T.y2.lasso.boot$t[,1],c(0.025,0.05,0.1,0.5,0.9,0.95,0.975),na.rm=TRUE);

save.image("./nhis13.RData");



#=======================================================================
library(stargazer);

lm.greg <- lm(update(model.GREG,"cancer~."),data=nhis13);
lm.obj <- lm(update(model,"cancer~."),data=nhis13);


sink("./models.tex");
stargazer(lm.greg,lm.obj,single.row=TRUE,report=c('vc*'),float=FALSE);
sink();






vars <- c("region","agegrp","gender","educ2","race","marst","employed","faminc_q");

props.samp1 <- NULL;
props.samp2 <- NULL;
props.pop <- NULL;
for (v in vars) {
  props.A <- (table(nhis13[,v])/nrow(nhis13));
  xtab.formula <- formula(sprintf("pwgtp~%s",v));
  props.B <- (xtabs(xtab.formula,data=acs13b)/sum(xtabs(xtab.formula,data=acs13b)));
  xtab.formula <- formula(sprintf("wtfa_sa~%s",v));
  props.C <- (xtabs(xtab.formula,data=nhis13)/sum(xtabs(xtab.formula,data=nhis13)));
 
 props.samp1 <- c(props.samp1,props.A);
 props.samp2 <- c(props.samp2,props.C);
 props.pop <- c(props.pop,props.B);
    cat(sprintf("%s %d %d\n",v,length(props.A),length(props.B)));

}

distr <- rbind(props.samp1,props.samp2,props.pop);
write.table(distr,file="./var.distr.csv",row.names=F,sep=",");

nhiswt<-c(nhis13$wts2)
gregwt2<-c(nhis13$wts3)
gregwt<-c(nhis13$wts4)
lassowt<-wts4

sqrt(var(nhiswt))
sqrt(var(gregwt))
sqrt(var(gregwt2))
sqrt(var(lassowt))

c(min(nhiswt),max(nhiswt))
c(min(gregwt),max(gregwt))
c(min(gregwt2),max(gregwt2))
c(min(lassowt),max(lassowt))

sqrt((21070498-19899327)^2+362833^2)
sqrt((20254449-19899327)^2+375064^2)
sqrt((20281603-19899327)^2+367900^2)
sqrt((20064671-19899327)^2+347586^2)


