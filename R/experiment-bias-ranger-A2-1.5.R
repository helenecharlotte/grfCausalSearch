### experiment-bias-ranger-A2-1.5.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jan 12 2021 (17:26) 
## Version: 
## Last-Updated: Jan 13 2021 (14:05) 
##           By: Thomas Alexander Gerds
##     Update #: 26
#----------------------------------------------------------------------
## 
### Commentary:
##
## Here we look into a bias of the IPCW.ranger 
## seen in a specific simulation setting.
##
## To illustrate this use the modified functions from github. We
## now fit one forest for A=0 and another for A=1 to calculate the weights.
##
## In addition num.tree=1 produces unbiased ATE results in this specific scenario.
## see below
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
## plot(f1,confint=FALSE,type="cuminc",xlim=c(0,1))
## abline(v=.5,col=2,lwd=2)

## crude observed differences -- these depend on effect.A2
f <- prodlim(Hist(time,delta)~A1+A2,data=d)
## plot(f,confint=FALSE,xlim=c(0,1))
## abline(v=.5,col=2,lwd=2)

#----------------------------------------------------------------------------------------------------
# 2. The remaining problem after running two forests for each intervention node:
#    one for A=0 and another for A=1 
#----------------------------------------------------------------------------------------------------
# for n=200

W <- runner(seed=179,cens=0.6,cores=25,
            effect.A2=1.5,M=1000,n=200,
            NT=150,method.weight="ranger",intervene="A2",
            formula.weight=Hist(time,delta)~A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,
            fit.separate=0,
            time.interest=.5,
            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49,min.node.size=NULL),
            verbose=1L)
## intervene   n    m time    cens  event1  event2 truth num.trees        mean        se    mean.se       bias   abs.bias coverage
##        A2 200 1000  0.5 0.10431 0.13466 0.40553     0       150 -0.02464767 0.0689244 0.07407326 0.02464767 0.05847687    0.956

# for n=1000
W <- runner(seed=179,cens=0.6,cores=25,
            effect.A2=1.5,M=1000,n=1000,
            NT=150,method.weight="ranger",intervene="A2",
            formula.weight=Hist(time,delta)~A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,
            fit.separate=0,
            time.interest=.5,
            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE,tuning.time=0.49,min.node.size=NULL),
            verbose=1L)
W
##   intervene    n    M time     cens   event1   event2 truth num.trees        mean         se    mean.se       bias   abs.bias coverage
##          A2 1000 1000  0.5 0.104546 0.135165 0.405256     0       150 -0.02471412 0.02978624 0.03244842 0.02471412 0.03135142    0.906

#----------------------------------------------------------------------------------------------------
# 3. IPCW.ranger with a single tree 
#----------------------------------------------------------------------------------------------------
Tree1 <- foreach(sample.size=c(200,500,750,1000),.combine="rbind")%dopar%{
    V <- runner(seed=179,cens=0.6,cores=25,
                effect.A2=1.5,M=1000,n=sample.size,
                NT=150,method.weight="ranger",intervene="A2",
                formula.weight=Hist(time,delta)~A1+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6,
                fit.separate=0,
                time.interest=0.5,
                args.weight=list(num.tree=1,replace=FALSE,probability=TRUE,tuning.time=0.49,min.node.size=NULL),
                verbose=1L)
    V
}
Tree1
# The bias gets small but seems to increase with the sample size:

##  intervene    n    M time      cens   event1   event2 truth num.trees         mean         se    mean.se        bias   abs.bias coverage
##         A2  200 1000  0.5 0.1043100 0.134660 0.405530     0       150 -0.008950055 0.10861225 0.10821788 0.008950055 0.08604933    0.969
##         A2  500 1000  0.5 0.1048840 0.134910 0.404350     0       150 -0.008258395 0.06885970 0.06969068 0.008258395 0.05453398    0.953
##         A2  750 1000  0.5 0.1046933 0.135024 0.404896     0       150 -0.009419005 0.05870082 0.05773295 0.009419005 0.04708312    0.950
##         A2 1000 1000  0.5 0.1045460 0.135165 0.405256     0       150 -0.010624006 0.04971568 0.05013392 0.010624006 0.04063588    0.958

######################################################################
### experiment-bias-ranger-A2-1.5.R ends here
