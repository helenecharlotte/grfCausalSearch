### discussion-bias-ranger.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 12 2021 (17:26) 
## Version: 
## Last-Updated: Jan 13 2021 (13:16) 
##           By: Thomas Alexander Gerds
##     Update #: 19
#----------------------------------------------------------------------
## 
### Commentary:
##
## Here we discuss a bias of the IPCW.ranger modus
## seen in a specific simulation setting.
## To illustrate this check in the modified functions
## and note that the results now contain the overall percentages
## of events 1 and 2 before the time intererst alongside
## the percentages right censored before the time interest.
## 
#----------------------------------------------------------------------
## 
### Code:
library(grid)
library(grf)
library(data.table)
library(gridExtra)
library(scales)
library(randomForestSRC)
library(ranger)
library(riskRegression)
library(prodlim)
library(survival)
library(parallel)
library(foreach)
library(doParallel)
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
try(setwd("/home/ifsv/jhl781/Dropbox/grfCausalSearch/"),silent=TRUE)
try(setwd("~/Dropbox/grfCausalSearch/"),silent=TRUE)
source("./R/sim-data.R")
source("./R/causalhunter.R")
source("./R/weighter.R")
#source("./R/weighter.new.R")
source("./R/runner.R")
library(data.table)
library(ggplot2)

#----------------------------------------------------------------------------------------------------
# 1. The setting
#----------------------------------------------------------------------------------------------------
effect.A2 <- 1.5
form.T1 <- function(X, A){
    -1.1+as.numeric(X[, 1])*0.2-as.numeric(X[, 3])*0.1-A[, 1]*1.5}
form.T2 <- function(X, A) {
    0.1-as.numeric(X[, 2])*0.4-as.numeric(X[, 1])*0.33+effect.A2*A[, 2]}
set.seed(1919)
d <- sim.data(10000,CR=2,form.T1=form.T1,form.T2=form.T2,shape.T1=0.8,shape.T2=0.8,C.shape=0.6,which.A=2)
d[,table(delta)/.N]
d[time<=.5,table(delta)/.N]
## counterfactual world world where cause 2 is eliminated (and no censoring occurs)
d[,T1:=T1.1*A2 + T1.0*(1-A2)]
d[,dummy:=1]
f1 <- prodlim(Hist(T1,dummy)~A1+A2,data=d)
## no effect of A2 on T1
plot(f1,confint=FALSE,type="cuminc",xlim=c(0,1))
abline(v=.5,col=2,lwd=2)

## crude observed differences -- these depend on effect.A2
f <- prodlim(Hist(time,delta)~A1+A2,data=d)
plot(f,confint=FALSE,xlim=c(0,1))
#plot(f,confint=FALSE,cause="2")
abline(v=.5,col=2,lwd=2)

#----------------------------------------------------------------------------------------------------
# 1. The problem
#----------------------------------------------------------------------------------------------------
W <- runner(seed=179,cens=0.6,cores=25,
            effect.A2=1.5,M=1000,n=200,
            NT=150,method.weight="ranger",intervene="A2",
            formula.weight=Hist(time,delta)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,
            fit.separate=0,
            time.interest=.5,
            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49,min.node.size=NULL),
            verbose=1L)
W

# There is a bias 
##    intervene time    cens  event1  event2 truth num.trees       mean         se    mean.se      bias   abs.bias coverage
##           A2  0.5 0.10565 0.13475 0.40405     0       150 -0.0518181 0.07288361 0.07533538 0.0518181 0.07033412    0.905

#----------------------------------------------------------------------------------------------------
# 2. Increasing the sample size
#----------------------------------------------------------------------------------------------------
# The bias gets smaller but is still substantial when the sample size increases
foreach(n=c(200,500,750,1000,1500))%dopar%{
    V <- runner(seed=179,cens=0.6,cores=25,
                effect.A2=1.5,M=1000,n=1000,
                NT=150,method.weight="ranger",intervene="A2",
                formula.weight=Hist(time,delta)~A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,
                fit.separate=0,
                time.interest=0.5,
                args.weight=list(num.tree=1,replace=FALSE,probability=TRUE,tuning.time=0.49,min.node.size=NULL),
                verbose=1L)
    V
}

## >    intervene time     cens   event1   event2 truth num.trees        mean         se    mean.se       bias   abs.bias
## 1:        A2  0.5 0.104546 0.135165 0.405256     0       150 -0.02471412 0.02978624 0.03244842 0.02471412 0.03135142
   ## coverage
## 1:    0.906

## intervene time   cens   event1  event2 truth num.trees        mean         se    mean.se       bias   abs.bias coverage
##        A2  0.5 0.1053 0.134195 0.40659     0       150 -0.03084236 0.03143825 0.03264133 0.03084236 0.03635274     0.86

#----------------------------------------------------------------------------------------------------
# 3. Increasing some more
#----------------------------------------------------------------------------------------------------
V <- runner(seed=179,cens=0.6,cores=15,
            effect.A2=1.5,M=200,n=1500,
            NT=150,method.weight="ranger",intervene="A2",
            formula.weight=Hist(time,delta)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,
            fit.separate=0,
            time.interest=0.5,
            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49,min.node.size=NULL),
            verbose=1L)
V

#----------------------------------------------------------------------------------------------------
# 4. Focus the mighty source of the problem
#----------------------------------------------------------------------------------------------------

set.seed(19)
d <- sim.data(1000,CR=2,form.T1=form.T1,form.T2=form.T2,shape.T1=0.1,shape.T2=1.8,C.shape=0.6,which.A=2)
f <- prodlim(Hist(time,delta)~A1+A2,data=d)
plot(f,confint=FALSE,xlim=c(0,1))
abline(v=.5,col=2,lwd=2)
d[,cens:=as.numeric(delta!=1)]


A <- causalhunter(formula=Hist(time,delta)~A1+intervene(A2)+X1+X2+X3+X4+X5+X6,
                  formula.weight=Hist(time,delta)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,
                  times=0.5,
                  data=d)

Ykm <- weighter(Hist(time,delta)~A1+A2,data=d,method="km",times=.5,num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49)
data.table("A2=1"=mean(Ykm[d$A2==1]),"A2=0"=mean(Ykm[d$A2==0]),ATE=mean(Ykm[d$A2==1])-mean(Ykm[d$A2==0]))

Y0 <- weighter(Hist(time,delta)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,data=d[A2==0],method="ranger",times=.5,num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49)
Y1 <- weighter(Hist(time,delta)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,data=d[A2==1],method="ranger",times=.5,num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49)
d[A2==0,Y:=Y0]
d[A2==1,Y:=Y1]
mean(Y0)-mean(Y1)
data.table("A2=1"=mean(Y[d$A2==1]),"A2=0"=mean(Y[d$A2==0]),ATE=mean(Y[d$A2==1])-mean(Y[d$A2==0]))


Y <- weighter(Hist(time,delta)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,data=d,method="ranger",times=.5,num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49)
data.table("A2=1"=mean(Y[d$A2==1]),"A2=0"=mean(Y[d$A2==0]),ATE=mean(Y[d$A2==1])-mean(Y[d$A2==0]))

Yr <- weighter(Hist(time,delta)~A1+A2,data=d,method="ranger",times=.5,num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49)
data.table("A2=1"=mean(Yr[d$A2==1]),"A2=0"=mean(Yr[d$A2==0]),ATE=mean(Yr[d$A2==1])-mean(Yr[d$A2==0]))

grf.A <- causal_forest(X=model.matrix(~-1+A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,model.frame(formula=~A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,data=d)),
                       Y=Y,
                       W=d$A2,num.tree=150)
average_treatment_effect(grf.A)

grf.A <- causal_forest(X=model.matrix(~-1+A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,model.frame(formula=~A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,data=d)),
                       Y=matrix(d[["Y"]],ncol=1),
                       W=d$A2,num.tree=350)
average_treatment_effect(grf.A)

grf.A <- causal_forest(X=model.matrix(~-1+A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,model.frame(formula=~A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,data=d)),
                       Y=Ykm,
                       W=d$A2,num.tree=150)
average_treatment_effect(grf.A)

plot(Y,Ykm)
plot(Yr,Ykm)
plot(d$Y,Ykm)
plot(d$Y,Y)

######################################################################
### discussion-bias-ranger.R ends here
