### simulation-study-for-article.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 27 2020 (10:20) 
## Version: 
## Last-Updated: Jun 29 2020 (08:47) 
##           By: Thomas Alexander Gerds
##     Update #: 27
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
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
source("~/research/SoftWare/grfCausalSearch/R/sim-data.R")
source("~/research/SoftWare/grfCausalSearch/R/causalhunter.R")
source("~/research/SoftWare/grfCausalSearch/R/weighter.R")
source("~/research/SoftWare/grfCausalSearch/R/runner.R")

if (system("echo $USER",inter=TRUE)=="tag")
    setwd("~/research/SoftWare/grfCausalSearch/simulation-results")

if (FALSE){
    results.net <- do.call("rbind",lapply(rev(c(200,300,400,500,600,700,800,900,1000)),function(sample.size){
        do.call("rbind",lapply(c(0.5,.7,1.5),function(cens){
            do.call("rbind",lapply(c(0,0.5,0.75,1,1.25,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,effect.A2=effect,M=1000,n=sample.size,NT=50,method.weight="ranger",args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),verbose=0L)
                cbind(n=sample.size,Effect.A2=effect,cens=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file="./results-net.rds")
}

results.net.weight.formula <- do.call("rbind",lapply(rev(c(200,300,400,500,600,700,800,900,1000)),function(sample.size){
    cat(sample.size,"\n")
    do.call("rbind",lapply(c(.5,1.25,1.5),function(effect){
        cat(effect,"\n")
        v <- runner(formula.weight=Hist(time,delta)~A2+X1+X2,seed=179,cens=1.5,cores=25,effect.A2=effect,M=1000,n=sample.size,NT=50,method.weight="ranger",args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),verbose=0L)
        cbind(n=sample.size,Effect.A2=effect,censpar=1.5,v)
    }))
}))

saveRDS(results.net.weight.formula,file="./results-net-weight-formula.rds")

######################################################################
### simulation-study-for-article.R ends here
