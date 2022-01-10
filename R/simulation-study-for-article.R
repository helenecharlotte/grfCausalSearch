### simulation-study-for-article.R --- 
#----------------------------------------------------------------------
## Author: Hely & Tag
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

#-- set working directory; 
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
try(setwd("/home/ifsv/jhl781/Dropbox/grfCausalSearch/"),silent=TRUE)
try(setwd("~/Dropbox/grfCausalSearch/"),silent=TRUE)

#-- source code; 
source("./R/sim-data.R")
source("./R/causalhunter.R")
source("./R/weighter.R")
#source("./R/weighter.new.R")
source("./R/runner.R")

#-- how many repetitions?
M <- 500#1000#500#499#500#99#501#502#101
results.net <- do.call("rbind",lapply(rev(c(1000)),function(sample.size){
    do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
        do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
            v <- runner(seed=179,cens=cens,cores=25,
                        effect.A2=effect,M=M,n=sample.size,
                        NT=200,method.weight="ranger",
                        args.weight=list(num.tree=1,replace=FALSE,probability=TRUE),
                        verbose=1L)
            cbind(Effect.A2=effect,cens.shape=cens,v)
        }))
    }))
}))
saveRDS(results.net,file=paste0("./simulation-results/results2-net-weight-correct-numtree1-M", M, ".rds"))


if (FALSE) {#test
    M <- 500
    runner(seed=288,cens=cens,cores=1,
           effect.A2=0.2,M=1,n=300,
           fit.separate=TRUE,
           #formula.weight=Hist(time,delta)~X2+X3+A2,
           NT=200,method.weight="ranger",
           args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
           effect="net",
           verbose=1L)
}

if (FALSE) {#test
    M <- 310#400
    cens <- 0.4
    test1 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    #formula.weight=Hist(time,delta)~X1+X2+A2,
                    fit.separate=FALSE,
                    NT=200,method.weight="ranger",
                    args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test1,file=paste0("./simulation-results/test1-new-not-separate-net-M", M, ".rds"))
}

if (FALSE) {#test
    M <- 250#400
    cens <- 0.4
    test1 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    #formula.weight=Hist(time,delta)~X1,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test1,file=paste0("./simulation-results/test1-net-M", M, ".rds"))
}

if (FALSE) {
    M <- 400
    cens <- 0.4
    test2 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    formula.weight=Hist(time,delta)~A2+X1+X2,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test2,file=paste0("./simulation-results/test2-net-M", M, ".rds"))
}


if (FALSE) {
    M <- 400
    cens <- 0.4
    test21 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    formula.weight=Hist(time,delta)~A2+X1+X2+X3,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test21,file=paste0("./simulation-results/test21-net-M", M, ".rds"))
}

if (FALSE) {
    M <- 400
    cens <- 0.4
    test22 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    formula.weight=Hist(time,delta)~A2+X1+X2+X3+X4,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test22,file=paste0("./simulation-results/test22-net-M", M, ".rds"))
}

if (FALSE) {
    M <- 400
    cens <- 0.4
    test23 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    formula.weight=Hist(time,delta)~A2+X1+X2+X3+X4+X5,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test23,file=paste0("./simulation-results/test23-net-M", M, ".rds"))
}



if (FALSE) {#test
    M <- 400
    cens <- 0.4
    test3 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    formula.weight=Hist(time,delta)~A2+X1+X2+X3+X4+X5+X6,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test3,file=paste0("./simulation-results/test3-net-M", M, ".rds"))
}

if (FALSE) {
    M <- 400
    cens <- 0.4
    test31 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    formula.weight=Hist(time,delta)~A2+X1+X2+X3+X4+X5+X6+A1,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test31,file=paste0("./simulation-results/test31-net-M", M, ".rds"))
}

if (FALSE) {
    M <- 400
    cens <- 0.4
    test32 <- runner(seed=179,cens=cens,cores=25,
                     effect.A2=0.5,M=M,n=1000,
                     formula.weight=Hist(time,delta)~A2+X2+X3+A1+A3+A4+A5+A6+A7+A9+A10,
                     NT=50,method.weight="ranger",
                     args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                     effect="net",
                     verbose=1L)
    saveRDS(test32,file=paste0("./simulation-results/test32-net-M", M, ".rds"))
}


if (FALSE) {
    M <- 400
    cens <- 0.4
    test4 <- runner(seed=179,cens=cens,cores=25,
                    effect.A2=0.5,M=M,n=1000,
                    formula.weight=Hist(time,delta)~A2+X1+X2+X3+X4+X5+X6+A1+A3+A4+A5,
                    NT=50,method.weight="ranger",
                    args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                    effect="net",
                    verbose=1L)
    saveRDS(test4,file=paste0("./simulation-results/test4-net-M", M, ".rds"))
}

#-- run simulations; 
if (FALSE){
    results.net <- do.call("rbind",lapply(rev(c(200,300,400,500,600,700,800,900,1000)),function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=50,method.weight="ranger",
                            args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file=paste0("./simulation-results/results1-net-M", M, ".rds"))
    saveRDS(results.net,file=paste0("./simulation-results/results1-net-weight-correct-M", M, ".rds"))
}

if (FALSE){
    results.net <- do.call("rbind",lapply(rev(c(1000)),function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="ranger",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file=paste0("./simulation-results/results3-net-M", M, ".rds"))
}


if (FALSE) {
    results.net <- do.call("rbind",lapply(rev(c(1000)),function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="ranger",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file=paste0("./simulation-results/results2-net-weight-correct-M", M, ".rds"))
}

if (TRUE) {
    results.net <- do.call("rbind",lapply(rev(c(1000)),function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0.2,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=15,
                            formula.weight=Hist(time,delta)~A2+X1+X2,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="km",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file=paste0("./simulation-results/results3-net-km-weight-correct-M", M, ".rds"))
}

if (FALSE) {
    results.net <- do.call("rbind",lapply(rev(c(1000)),function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            formula.weight=Hist(time,delta)~A2+X1+X2,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="ranger",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file=paste0("./simulation-results/results3-net-weight-correct-formula-M", M, ".rds"))
}

if (FALSE) {
    results.net.weight.correct <- do.call("rbind",lapply(1000,function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="ranger",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            effect="crude",
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net.weight.correct,file=paste0("./simulation-results/results2-crude-weight-correct-M", M, ".rds"))
}

if (FALSE) {
    results.net.weight.correct <- do.call("rbind",lapply(1000,function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            formula.weight=Hist(time,delta)~A2+X1+X2,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="ranger",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            effect="crude",
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net.weight.correct,file=paste0("./simulation-results/results3-crude-weight-correct-formula-M", M, ".rds"))
}


if (FALSE) {
    results.net.weight.misspecified <- do.call("rbind",lapply(1000,function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(formula.weight=Hist(time,delta)~X1,
                            seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="ranger",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net.weight.misspecified,file=paste0("./simulation-results/results3-net-weight-misspecified-M", M, ".rds"))
}

if (FALSE) {
    results.net.weight.misspecified.A2 <- do.call("rbind",lapply(1000,function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(formula.weight=Hist(time,delta)~A2,
                            seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="ranger",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net.weight.misspecified.A2,
            file=paste0("./simulation-results/results3-net-weight-misspecified-A2-M", M, ".rds"))
}


if (FALSE) {
    results.crude.weight.misspecified <- do.call("rbind",lapply(1000,function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(formula.weight=Hist(time,delta)~X1,
                            seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="ranger",
                            effect="crude",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.crude.weight.misspecified,file=paste0("./simulation-results/results3-crude-weight-misspecified-M", M, ".rds"))
}

if (FALSE) {
    results.crude.weight.misspecified.A2 <- do.call("rbind",lapply(1000,function(sample.size){
        do.call("rbind",lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind",lapply(c(-0.2,0,0.2,0.5,0.7,1.5),function(effect){
                v <- runner(formula.weight=Hist(time,delta)~A2,
                            seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=150,method.weight="ranger",
                            effect="crude",
                            args.weight=list(num.tree=150,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.crude.weight.misspecified.A2,
            file=paste0("./simulation-results/results3-crude-weight-misspecified-A2-M", M, ".rds"))
}


######################################################################
### simulation-study-for-article.R ends here
