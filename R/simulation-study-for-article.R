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
library(nleqslv)

#-- set working directory; 
if (system("echo $USER",inter=TRUE)=="tag") {
    setwd("~/research/SoftWare/grfCausalSearch/")
} else if (system("echo $USER",intern=TRUE)%in%c("jhl781")) {
    setwd("/home/ifsv/jhl781/Dropbox/PhD/grfCausalSearch/")
} else if (system("echo $USER",intern=TRUE)%in%c("helene")) {
    setwd("~/Dropbox/PhD/grfCausalSearch/")
}

#-- source code; 
source("./R/sim-data.R")
source("./R/causalhunter.R")
source("./R/weighter.R")
#source("./R/weighter.new.R")
source("./R/runner.R")

#-- how many repetitions?
M <- 500

#-------------- FIGURE 1 -------------#

if (TRUE) {
    results.net.a <- do.call("rbind", lapply(rev(c(1000)),function(sample.size){
        do.call("rbind", lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~X1+X2+A2,
                            effect="net",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net.a,file=paste0("./simulation-results/results-net-km-weight-a-M", M, ".rds"))
    results.net.b <- do.call("rbind", lapply(rev(c(1000)),function(sample.size){
        do.call("rbind", lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~A2,
                            effect="net",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net.b,file=paste0("./simulation-results/results-net-km-weight-b-M", M, ".rds"))
    results.net.c <- do.call("rbind", lapply(rev(c(1000)),function(sample.size){
        do.call("rbind", lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~1,
                            effect="net",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net.c,file=paste0("./simulation-results/results-net-km-weight-c-M", M, ".rds"))
}


if (TRUE) {
    results.crude.a <- do.call("rbind", lapply(rev(c(1000)),function(sample.size){
        do.call("rbind", lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~X1+X2+A2,
                            effect="crude",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.crude.a,file=paste0("./simulation-results/results-crude-km-weight-a-M", M, ".rds"))
    results.crude.b <- do.call("rbind", lapply(rev(c(1000)),function(sample.size){
        do.call("rbind", lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~A2,
                            effect="crude",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.crude.b,file=paste0("./simulation-results/results-crude-km-weight-b-M", M, ".rds"))
    results.crude.c <- do.call("rbind", lapply(rev(c(1000)),function(sample.size){
        do.call("rbind", lapply(c(0.2,0.4,0.6),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~1,
                            effect="crude",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.crude.c,file=paste0("./simulation-results/results-crude-km-weight-c-M", M, ".rds"))
}

#-------------- FIGURE 2 -------------#

if (TRUE) {
    results.net <- do.call("rbind", lapply(rev(c(200,300,400,500,750,1000,1500,2000)),function(sample.size){
        do.call("rbind", lapply(c(0.4),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~X1+X2+A2,
                            effect="net",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file=paste0("./simulation-results/results-hely-net-km-weight-correct-M", M, ".rds"))
    results.net <- do.call("rbind", lapply(rev(c(200,300,400,500,750,1000,1500,2000)),function(sample.size){
        do.call("rbind", lapply(c(0.4),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~1,
                            effect="net",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.net,file=paste0("./simulation-results/results-hely-net-km-weight-unadjusted-M", M, ".rds"))
}

if (TRUE) {
    results.crude <- do.call("rbind", lapply(rev(c(200,300,400,500,750,1000,1500,2000)),function(sample.size){
        do.call("rbind", lapply(c(0.4),function(cens){
            do.call("rbind", lapply(c(0.5, 1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~X1+X2+A2,
                            effect="crude",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.crude,file=paste0("./simulation-results/results-hely-crude-km-weight-correct-M", M, ".rds"))
    results.crude <- do.call("rbind", lapply(rev(c(200,300,400,500,750,1000,1500,2000)),function(sample.size){
        do.call("rbind", lapply(c(0.4),function(cens){
            do.call("rbind", lapply(c(-0.5,0.5,1.5),function(effect){
                v <- runner(seed=179,cens=cens,cores=25,
                            effect.A2=effect,M=M,n=sample.size,
                            NT=200,method.weight="km",
                            formula.weight=Hist(time,delta)~1,
                            effect="crude",
                            args.weight=list(num.tree=200,replace=FALSE,probability=TRUE),
                            verbose=1L)
                cbind(n=sample.size,Effect.A2=effect,cens.shape=cens,v)
            }))
        }))
    }))
    saveRDS(results.crude,file=paste0("./simulation-results/results-hely-crude-km-weight-unadjusted-M", M, ".rds"))
}


######################################################################
### simulation-study-for-article.R ends here
