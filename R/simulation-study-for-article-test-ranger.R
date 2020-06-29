#-------------------------------------------------------------------------------------------#
## set working directory
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){
    setwd("/home/ifsv/jhl781/research/phd/random-forest/r-code/causal-search-dst/")
} else {
    setwd("~/research/phd/random-forest/r-code/causal-search-dst/")
}

library(grid)
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
 
source("master.R")

source("plot-output-fun.R")
source("grfcens.R")
source("compute-aalen-johansen.R")

#--------------------- test with simulations ------------------------#'

n <- 200#1000
M <- 150
n.A <- 10

censoring <- "medium" #"medium" # "max", "min"

compute.truth <- FALSE

parameter <- "net" 

if (censoring=="medium") {
    C.shape <- 0.4
} else if (censoring=="max") {
    C.shape <- 0.2
} else if (censoring=="min") {
    C.shape <- 0.6
}



#-- compute true values for all treatments and both theta1 and theta2:

if (compute.truth) {
    psi0.list.V.1 <- unlist(lapply(1:n.A, function(a) {
        sim.data(1e6, compute.psi=1, CR=2, which.A=a,
                 form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                          as.numeric(X[, 1])*0.33 +
                                          1.5*A[, 2])}))
    names(psi0.list.V.1) <- paste0("A", 1:10)
    saveRDS(psi0.list.V.1, file=paste0("output/psi0.list.V.1", ".rds"))

    psi0.list.V.2 <- unlist(lapply(1:n.A, function(a) {
        sim.data(1e6, compute.psi=2, CR=2, which.A=a,
                 form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                          as.numeric(X[, 1])*0.33 +
                                          1.5*A[, 2])}))
    names(psi0.list.V.2) <- paste0("A", 1:10)
    saveRDS(psi0.list.V.2, file=paste0("output/psi0.list.V.2", ".rds"))
}

#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 15
} else {
    no_cores <- detectCores() - 1
}


registerDoParallel(no_cores)

ranger.list.2 <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
                         ) %dopar% {
                             print(m)
                             dt <- sim.data(n, CR=2, seed=201010120+m,
                                            C.shape=C.shape,
                                            form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                                                     as.numeric(X[, 1])*0.33 +
                                                                     1.5*A[, 2])
                             if (parameter=="crude") {
                                 return(do.call("rbind", lapply(1:10, function(a)
                                     grfcens(as.formula(paste0("Hist(time,delta)~intervene(A",a,")+X1+X2+X3+X4+X5+X6")),
                                             method.weight="ranger",
                                             data=dt, times=0.5))))
                             } else {
                                 return(do.call("rbind", lapply(1:10, function(a)
                                     grfcens(as.formula(paste0("Hist(time,delta)~intervene(A",a,")+X1+X2+X3+X4+X5+X6")),
                                             method.weight="ranger",
                                             CR.as.censoring=TRUE, 
                                             data=dt, times=0.5))))
                             }
                         }

stopImplicitCluster()

if (parameter=="crude") {
    saveRDS(ranger.list.2, file=paste0("output/ATE-list-fast",
                                       "-V-theta2-ranger-fast-n",n,"-M",M,
                                       "-censoring", censoring, ".rds"))
} else {
    saveRDS(ranger.list.2, file=paste0("output/ATE-list-fast",
                                       "-V-theta1-ranger-fast-n",n,"-M",M,
                                       "-censoring", censoring, ".rds"))
}

if (FALSE) {

    n <- 400; M <- 1000; n.A <- 10
    #M <- 1000

    new <- "-censoring"
    censoring <- "medium"
    parameter <- "net"

    if (parameter=="crude") {
        length(ranger.list.2 <- readRDS(file=paste0("output/ATE-list-fast",
                                                    "-V-theta2-ranger-fast-n",n,"-M",M,"-censoring",
                                                    censoring, ".rds")))
        psi0.list.V.2 <- readRDS(file=paste0("output/psi0.list.V.2", ".rds"))

        cbind(mean.fun2(ranger.list.2), ranger.list.2[[1]])
        coverage.fun2(ranger.list.2)
    } else {
        length(ranger.list.1 <- readRDS(file=paste0("output/ATE-list-fast",
                                                    "-V-theta1-ranger-fast-n",n,"-M",M,"-censoring",
                                                    censoring, ".rds")))
        psi0.list.V.1 <- readRDS(file=paste0("output/psi0.list.V.1", ".rds"))

        print(cbind(mean.fun2(ranger.list.1), ranger.list.1[[1]], psi0=psi0.list.V.1))
        print(coverage.fun2(ranger.list.1, psi0=psi0.list.V.1))
    }
    
}

coverage.fun2 <- function(ATE.list, psi0=psi0.list.V.2) {
    return(coverage.1 <- unlist(lapply(1:n.A, function(a) {
        return(mean(unlist(lapply(ATE.list, function(est) {
            #est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
            est.a <- est[a, ]
            return(as.numeric(na.omit(est.a["lower"]<=psi0[a] & psi0[a]<=est.a["upper"])))
        }))))
    })))
}

mean.fun2 <- function(ATE.list, psi0=psi0.list.V.2) {
    return(coverage.1 <- unlist(lapply(1:n.A, function(a) {
        return(mean(unlist(lapply(ATE.list, function(est) {
            #est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
            est.a <- est[a, "ate"]
            return(as.numeric(est.a))
        }))))
    })))
}


if (FALSE) {
for (m in 1:M) {

    #---- scenario V:
    dt <- sim.data(n, CR=2, seed=201010120+m,
                   C.shape=C.shape,
                   form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                            as.numeric(X[, 1])*0.33 +
                                            1.5*A[, 2])

    if (FALSE) {

        nevents[[m]] <- table(dt[time<=0.5, delta])

        saveRDS(nevents, file=paste0("output/nevents",
                                     "-V-n",n,"-M",M,"-censoring", censoring, ".rds"))
        
        # estimate theta1:
        ATE.list.V.1.unadj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE, 
                                                   X.vars=paste0("X", 1:10))
        ATE.list.V.1.adjA[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                                  X.vars=paste0("X", 1:10),#method.weight="ranger", 
                                                  stratify.CR="A2")
        ATE.list.V.1.adjA2[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                                   X.vars=paste0("X", 1:10),
                                                   stratify.A=TRUE)
        ATE.list.V.1.adj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                                 X.vars=paste0("X", 1:10),
                                                 stratify.CR="A2+X1+X2")
   
        # estimate theta2:
        ATE.list.V.2.unadj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                                   X.vars=paste0("X", 1:10))
        ATE.list.V.2.adjA[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                                  X.vars=paste0("X", 1:10),
                                                  stratify.CR="A2")
        ATE.list.V.2.adjA2[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                                   X.vars=paste0("X", 1:10),
                                                   stratify.A=TRUE)
        ATE.list.V.2.adj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                                 X.vars=paste0("X", 1:10),
                                                 stratify.CR="A2+X1+X2")

        #-- save output:
    
        saveRDS(ATE.list.V.1.unadj, file=paste0("output/ATE-list-fast",
                                                "-V-theta1-unadj-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
        saveRDS(ATE.list.V.1.adjA, file=paste0("output/ATE-list-fast",
                                               "-V-theta1-adjA-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
        saveRDS(ATE.list.V.1.adjA2, file=paste0("output/ATE-list-fast",
                                                "-V-theta1-adjA2-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
        saveRDS(ATE.list.V.1.adj, file=paste0("output/ATE-list-fast",
                                              "-V-theta1-adj-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
        saveRDS(ATE.list.V.2.unadj, file=paste0("output/ATE-list-fast",
                                                "-V-theta2-unadj-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
        saveRDS(ATE.list.V.2.adjA, file=paste0("output/ATE-list-fast",
                                               "-V-theta2-adjA-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
        saveRDS(ATE.list.V.2.adjA2, file=paste0("output/ATE-list-fast",
                                                "-V-theta2-adjA2-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
        saveRDS(ATE.list.V.2.adj, file=paste0("output/ATE-list-fast",
                                              "-V-theta2-adj-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))

    }
    
    #ranger.list.1 <- do.call("rbind", lapply(1:10, function(a)
    #    grfcens(as.formula(paste0("Hist(time,delta==2)~intervene(A",a,")+X1+X2+X3+X4+X5+X6")),
    #            data=dt,times=10)))
    ranger.list.2[[m+1]] <- do.call("rbind", lapply(1:10, function(a)
        grfcens(as.formula(paste0("Hist(time,delta)~intervene(A",a,")+X1+X2+X3+X4+X5+X6")),
                method.weight="ranger",
                data=dt, times=0.5)))
    saveRDS(ranger.list.2, file=paste0("output/ATE-list-fast",
                                       "-V-theta2-ranger-fast-n",n,"-M",M,"-censoring", censoring, ".rds"))
    #-- test with forest:

    if (FALSE) {
        rsf.LR <- rfsrc(formula(paste0("Surv(time, delta) ~ ", paste0("X", 1:6, collapse="+"), "+",
                                       paste0("A", 1:10, collapse="+"))),
                        importance=TRUE,
                        cause=1,
                        splitrule = "logrank",
                        data=dt)

        rsf.LR.ranked[[m]] <- sort(100 * rsf.LR$importance[ ,1])

        saveRDS(rsf.LR.ranked, file=paste0("output/ATE-list-fast",
                                           "-V-RSF-LR-ranked-n",n,"-M",M,"-censoring", censoring, ".rds"))
    
        rsf.gray <- rfsrc(formula(paste0("Surv(time, delta) ~ ", paste0("X", 1:6, collapse="+"), "+",
                                         paste0("A", 1:10, collapse="+"))),
                          importance=TRUE,
                          cause=1,
                          splitrule = "logrankCR",
                          data=dt)

        rsf.gray.ranked[[m]] <- sort(100 * rsf.gray$importance[ ,1])
    
        saveRDS(rsf.gray.ranked, file=paste0("output/ATE-list-fast",
                                             "-V-RSF-gray-ranked-n",n,"-M",M,"-censoring", censoring, ".rds"))
    }
    
}
}

print("done here")

coverage.fun <- function(ATE.list) {
    psi0 <- ATE.list[[1]]
    return(coverage.1 <- unlist(lapply(1:n.A, function(a) {
        return(mean(unlist(lapply(ATE.list[-1], function(est) {
            est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
            est.a <- est[a, ]
            return(as.numeric(na.omit(est.a$CI.lwr<=psi0[a] & psi0[a]<=est.a$CI.upr)))
        }))))
    })))
}

mean.fun <- function(ATE.list) {
    psi0 <- ATE.list[[1]]
    return(coverage.1 <- unlist(lapply(1:n.A, function(a) {
        return(mean(unlist(lapply(ATE.list[-1], function(est) {
            est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
            est.a <- est[a, ]
            return(as.numeric(est.a))
        }))))
    })))
}



if (FALSE) {

    n <- 1000; M <- 1000; n.A <- 10
    #M <- 1000

    new <- "-censoring"
    censoring <- "medium"


    length(rsf.LR.ranked <- readRDS(file=paste0("output/ATE-list-fast",
                                                "-V-RSF-LR-ranked-n",n,"-M",M,
                                                "-censoring", censoring, ".rds")))

    length(rsf.gray.ranked <- readRDS(file=paste0("output/ATE-list-fast",
                                                  "-V-RSF-gray-ranked-n",n,"-M",M,
                                                  "-censoring", censoring, ".rds")))

    A1.ranked.LR <- unlist(lapply(rsf.LR.ranked, function(out) {
        (1:length(out))[names(out[rev(order(out))])=="A1"]
    }))

    A2.ranked.LR <- unlist(lapply(rsf.LR.ranked, function(out) {
        (1:length(out))[names(out[rev(order(out))])=="A2"]
    }))

    A1.ranked.gray <- unlist(lapply(rsf.gray.ranked, function(out) {
        (1:length(out))[names(out[rev(order(out))])=="A1"]
    }))

    A2.ranked.gray <- unlist(lapply(rsf.gray.ranked, function(out) {
        (1:length(out))[names(out[rev(order(out))])=="A2"]
    }))

    ranked.rsf <- data.frame(A1.ranked.LR=A1.ranked.LR,
                             A1.ranked.gray=A1.ranked.gray,
                             A2.ranked.LR=A2.ranked.LR,
                             A2.ranked.gray=A2.ranked.gray)

    mean(ranked.rsf$A1.ranked.LR==1)
    mean(ranked.rsf$A1.ranked.gray==1)
    sd(ranked.rsf$A1.ranked.LR)
    sd(ranked.rsf$A1.ranked.gray)
    mean(ranked.rsf$A2.ranked.LR==1)
    mean(ranked.rsf$A2.ranked.gray==1)

    par(mfrow=c(2,1))
    hist(ranked.rsf$A2.ranked.LR)
    hist(ranked.rsf$A2.ranked.gray)
    par(mfrow=c(1,1))    
}

if (FALSE) {


    n <- 1000; M <- 100; n.A <- 10
    #M <- 1000

    new <- "-censoring"
    censoring <- "medium"

    length(ranger.list.2 <- readRDS(file=paste0("output/ATE-list-fast",
                                                "-V-theta2-ranger-fast-n",n,"-M",M,"-censoring",
                                                censoring, ".rds")))

    cbind(mean.fun2(ranger.list.2), ranger.list.2[[1]])
    coverage.fun2(ranger.list.2)
}

if (FALSE) {

    n <- 1000; M <- 1000; n.A <- 10
    #M <- 1000

    new <- "-censoring"
    censoring <- "medium"

    #--- get output:

    length(nevents <- readRDS(file=paste0("output/nevents",
                                          "-V-n",n,"-M",M,"-censoring", censoring, ".rds")))

    mean(unlist(lapply(nevents, function(nevent) {nevent[1]}))) # no of censoring events
    mean(unlist(lapply(nevents, function(nevent) {nevent[2]}))) # no of events of interest
    mean(unlist(lapply(nevents, function(nevent) {nevent[3]}))) # no of competing events
    
    length(ATE.list.V.1.unadj <- readRDS(file=paste0("output/ATE-list-fast",
                                                     "-V-theta1-unadj-fast-n",n,"-M",M,new, censoring, ".rds")))
    length(ATE.list.V.1.adj <- readRDS(file=paste0("output/ATE-list-fast",
                                                   "-V-theta1-adj-fast-n",n,"-M",M,new, censoring, ".rds")))
    length(ATE.list.V.1.adjA <- readRDS(file=paste0("output/ATE-list-fast",
                                                    "-V-theta1-adjA-fast-n",n,"-M",M,new, censoring, ".rds")))
    length(ATE.list.V.1.adjA2 <- readRDS(file=paste0("output/ATE-list-fast",
                                                     "-V-theta1-adjA2-fast-n",n,"-M",M,new, censoring, ".rds")))
    length(ATE.list.V.2.unadj <- readRDS(file=paste0("output/ATE-list-fast",
                                                     "-V-theta2-unadj-fast-n",n,"-M",M,new, censoring, ".rds")))
    length(ATE.list.V.2.adj <- readRDS(file=paste0("output/ATE-list-fast",
                                                   "-V-theta2-adj-fast-n",n,"-M",M,new, censoring, ".rds")))
    length(ATE.list.V.2.adjA <- readRDS(file=paste0("output/ATE-list-fast",
                                                    "-V-theta2-adjA-fast-n",n,"-M",M,new, censoring, ".rds")))
    length(ATE.list.V.2.adjA2 <- readRDS(file=paste0("output/ATE-list-fast",
                                                    "-V-theta2-adjA2-fast-n",n,"-M",M,new, censoring, ".rds")))    
 
    coverage.fun(ATE.list.V.1.unadj)
    coverage.fun(ATE.list.V.1.adj)
    coverage.fun(ATE.list.V.1.adjA)
    coverage.fun(ATE.list.V.1.adjA2)

    coverage.fun(ATE.list.V.2.unadj)
    coverage.fun(ATE.list.V.2.adj)
    coverage.fun(ATE.list.V.2.adjA)
    coverage.fun(ATE.list.V.2.adjA2)

    cbind(mean.fun(ATE.list.V.2.adjA), ATE.list.V.2.adjA[[1]])
    
    ATE.list.V.1.unadj[[1]]
    ATE.list.V.2.unadj[[1]]

    mean(unlist(lapply(2:length(ATE.list.V.1.adjA),
                       function(jj) {ATE.list.V.1.adjA[[jj]]["A2", "estimate"]})))
    mean(unlist(lapply(2:length(ATE.list.V.1.adjA2),
                       function(jj) {ATE.list.V.1.adjA2[[jj]]["A2", "estimate"]})))                     
    mean(unlist(lapply(2:length(ATE.list.V.1.unadj),
                       function(jj) {ATE.list.V.1.unadj[[jj]]["A2", "estimate"]})))

    mean(unlist(lapply(2:length(ATE.list.V.2.adjA),
                       function(jj) {ATE.list.V.2.adjA[[jj]]["A2", "estimate"]})))
    mean(unlist(lapply(2:length(ATE.list.V.2.adjA2),
                       function(jj) {ATE.list.V.2.adjA2[[jj]]["A2", "estimate"]})))                     
    mean(unlist(lapply(2:length(ATE.list.V.2.unadj),
                       function(jj) {ATE.list.V.2.unadj[[jj]]["A2", "estimate"]})))
    
    source("plot-output-fun.R")

    margin <- theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
   
    (p.V.1.unadj <- plot.output(ATE.list.V.1.unadj, n.A=3, pos3=0.25, 
                                sim=TRUE, text.pos=0.019, psi=1, text.pos.cov=0.02,
                                which="net",
                                title=bquote(hat(bar(theta))[net] ~ "(unadjusted weights)")))

    (p.V.2.unadj <- plot.output(ATE.list.V.2.unadj, n.A=3, pos3=0.25, which="crude", 
                                sim=TRUE, text.pos=0.016, psi=1:2, text.pos.cov=0.016,
                                title=bquote(hat(bar(theta))[crude] ~ "(unadjusted weights)")))

    p.V.12.unadj <- arrangeGrob(grobs=lapply(list(p.V.1.unadj, p.V.2.unadj), "+", margin),
                                nrow=1)

    ggsave(paste0("~/research/phd/random-forest/worg/figures/", "p-V-12-unadj-n", n, "-M", M,new, censoring, ".pdf"),
           p.V.12.unadj,
           width=5*2, height=3.5)

    ggsave(paste0("~/Dropbox/PhD/org-helene/figures/", "p-V-12-unadj-n", n, "-M", M,new, censoring, ".pdf"),
           p.V.12.unadj,
           width=5*2, height=3.5)

        
    
    (p.V.1.adjA <- plot.output(ATE.list.V.1.adjA, n.A=3, pos3=0.25, 
                               sim=TRUE, text.pos=0.019, psi=1, text.pos.cov=0.02,
                               which="net",
                               title=bquote(hat(bar(theta))[net] ~ "(weights adjusted for" ~ A[2]*")")))

    (p.V.2.adjA <- plot.output(ATE.list.V.2.adjA, n.A=3, pos3=0.25, which="crude",
                               sim=TRUE, text.pos=0.016, psi=1:2, text.pos.cov=0.016,
                               title=bquote(hat(bar(theta))[crude] ~ "(weights adjusted for" ~ A[2]*")")))

    p.V.12.adjA <- arrangeGrob(grobs=lapply(list(p.V.1.adjA, p.V.2.adjA), "+", margin),
                               nrow=1)

    ggsave(paste0("~/research/phd/random-forest/worg/figures/", "p-V-12-adjA-n", n, "-M", M,new, censoring, ".pdf"),
           p.V.12.adjA,
           width=5*2, height=3.5)

    ggsave(paste0("~/Dropbox/PhD/org-helene/figures/", "p-V-12-adjA-n", n, "-M", M,new, censoring, ".pdf"),
           p.V.12.adjA,
           width=5*2, height=3.5)

       
    (p.V.1.adj <- plot.output(ATE.list.V.1.adj, n.A=3, pos3=0.25,
                               sim=TRUE, text.pos=0.019, psi=1, text.pos.cov=0.02,
                               which="net",
                               title=bquote(hat(bar(theta))[net] ~ "(weights adjusted for" ~ A[2]*","*X[1]*","*X[2]*")")))

    (p.V.2.adj <- plot.output(ATE.list.V.2.adj, n.A=3, pos3=0.25,
                              which="crude",
                              sim=TRUE, text.pos=0.016, psi=1:2, text.pos.cov=0.016,
                              title=bquote(hat(bar(theta))[crude] ~ "(weights adjusted for" ~ A[2]*","*X[1]*","*X[2]*")")))

    p.V.12.adj <- arrangeGrob(grobs=lapply(list(p.V.1.adj, p.V.2.adj), "+", margin),
                              nrow=1)

    ggsave(paste0("~/research/phd/random-forest/worg/figures/", "p-V-12-adj-n", n, "-M", M,new, censoring, ".pdf"),
           p.V.12.adj,
           width=5*2, height=3.5)

    ggsave(paste0("~/Dropbox/PhD/org-helene/figures/", "p-V-12-adj-n", n, "-M", M,new, censoring, ".pdf"),
           p.V.12.adj,
           width=5*2, height=3.5)


     


}


if (FALSE) { ### how often is A1 ranked nr 1?

    read.fun <- function(which="V-theta1-adj", n=1000, M=500, which.A=1,
                         new="-censoring",
                         censoring="medium") {
        ATE.list <- readRDS(file=paste0("output/ATE-list-fast",
                                        "-", which, "-fast-n",n,"-M",M,new, censoring, ".rds"))
        ranked.1 <- unlist(lapply(ATE.list[-1], function(ate) {
            rownames(ate)[order(ate[, "estimate"])][1]
        }))
        return(ranked.1)
    }

    #read.fun()

    ranked.theta1.A1 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(n=sample.size)=="A1")
        }))

    ranked.theta2.A1 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adj", n=sample.size)=="A1")
        }))

    ranked.theta1.A2 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(n=sample.size)=="A2")
        }))

    ranked.theta2.A2 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adj", n=sample.size)=="A2")
        }))

    ranked.theta1.A3 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(n=sample.size)=="A3")
        }))

    ranked.theta2.A3 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adj", n=sample.size)=="A3")
        }))

     ranked.theta1.A4 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(n=sample.size)=="A4")
        }))

    ranked.theta2.A4 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adj", n=sample.size)=="A4")
        }))

     ranked.theta1.A5 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(n=sample.size)=="A5")
        }))

    ranked.theta2.A5 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adj", n=sample.size)=="A5")
        }))

     ranked.theta1.A6 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(n=sample.size)=="A6")
        }))

    ranked.theta2.A6 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adj", n=sample.size)=="A6")
        }))

    ranked.theta1.A7 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(n=sample.size)=="A7")
        }))

    ranked.theta2.A7 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adj", n=sample.size)=="A7")
        }))

    ## TEST FOR DEFENSE: 
    sum(unlist(lapply(1:10, function(a) {
        mean(read.fun(which="V-theta2-adj", n=100)==paste0("A",a))
    })))

    sum(unlist(lapply(1:10, function(a) {
        mean(read.fun(which="V-theta1-adj", n=100)==paste0("A",a))
    })))

    ranked.all.1 <- data.frame(do.call("rbind",
                                       lapply(c(100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                           unlist(lapply(1:10, function(a) {
                                               mean(read.fun(which="V-theta1-adj", n=sample.size)==paste0("A",a))
                                           }))
                                       })))

    names(ranked.all.1) <- paste0("A", 1:10)
    ranked.all.1$sample.size <- c(100, 200, 500, 1000, 1500, 2000)

    ranked.all.2 <- data.frame(do.call("rbind",
                                       lapply(c(100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                           unlist(lapply(1:10, function(a) {
                                               mean(read.fun(which="V-theta2-adj", n=sample.size)==paste0("A",a))
                                           }))
                                       })))

    names(ranked.all.2) <- paste0("A", 1:10)
    ranked.all.2$sample.size <- c(100, 200, 500, 1000, 1500, 2000) 

    ranked.all.11 <- melt(ranked.all.1, id="sample.size")
    ranked.all.22 <- melt(ranked.all.2, id="sample.size")

    ranked.all.11$parameter <- "1"
    ranked.all.22$parameter <- "2"
    
    ranked.all <- rbind(ranked.all.11, ranked.all.22)
    
    (p.all <- ggplot(ranked.all) +
         geom_line(aes(x=sample.size, y=value, col=variable)) + theme_bw() +
         theme(legend.title=element_blank()) +
         facet_wrap(parameter ~ .) +
         xlab("sample size") + ylab("fraction of times ranked no. 1"))

    ggsave(paste0("~/research/phd/random-forest/worg/figures/",
                  "plot-for-defense", ".pdf"),
           p.all,
           width=6.5, height=5)
    
    ranked.theta <- data.frame(censoring="adjusted weights (a)",#"correctly~adjusted~weights",
                               sample.size=c(100, 200, 500, 1000, 1500, 2000),
                               treatment=c(rep("A[1]~(with~direct~effect~on~event~of~interest)", length(ranked.theta1.A1)),
                                           rep("A[2]~(with~effect~only~on~the~competing~risk~event)", length(ranked.theta1.A1)),
                                           rep("A[3]~(with~no~effect)", length(ranked.theta1.A1))),
                               theta=c(rep("bar(theta)[net]", 3*length(ranked.theta1.A1)),
                                       rep("bar(theta)[crude]", 3*length(ranked.theta2.A1))),
                               rank=c(ranked.theta1.A1,
                                      ranked.theta1.A2,
                                      ranked.theta1.A3,
                                      ranked.theta2.A1,
                                      ranked.theta2.A2,
                                      ranked.theta2.A3))

    ranked.theta1.A1 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta1-adjA", n=sample.size)=="A1")
        }))

    ranked.theta2.A1 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adjA", n=sample.size)=="A1")
        }))

     ranked.theta1.A2 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta1-adjA", n=sample.size)=="A2")
        }))

    ranked.theta2.A2 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adjA", n=sample.size)=="A2")
        }))

     ranked.theta1.A3 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta1-adjA", n=sample.size)=="A3")
        }))

    ranked.theta2.A3 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-adjA", n=sample.size)=="A3")
        }))
    
    ranked.theta.censoring.max <- data.frame(censoring="adjusted weights (b)",
                                             sample.size=c(100, 200, 500, 1000, 1500, 2000),
                                             treatment=c(rep("A[1]~(with~direct~effect~on~event~of~interest)", length(ranked.theta1.A1)),
                                                         rep("A[2]~(with~effect~only~on~the~competing~risk~event)", length(ranked.theta1.A1)),
                                                         rep("A[3]~(with~no~effect)", length(ranked.theta1.A1))),
                                             theta=c(rep("bar(theta)[net]", 3*length(ranked.theta1.A1)),
                                                     rep("bar(theta)[crude]", 3*length(ranked.theta2.A1))),
                                             rank=c(ranked.theta1.A1,
                                                    ranked.theta1.A2,
                                                    ranked.theta1.A3,
                                                    ranked.theta2.A1,
                                                    ranked.theta2.A2,
                                                    ranked.theta2.A3))


    ranked.theta1.A1 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta1-unadj", n=sample.size)=="A1")
        }))

    ranked.theta2.A1 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-unadj", n=sample.size)=="A1")
        }))

     ranked.theta1.A2 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta1-unadj", n=sample.size)=="A2")
        }))

    ranked.theta2.A2 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-unadj", n=sample.size)=="A2")
        }))

     ranked.theta1.A3 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta1-unadj", n=sample.size)=="A3")
        }))

    ranked.theta2.A3 <- unlist(lapply(c(#100,
        100, 200, 500, 1000, 1500, 2000), function(sample.size) {
            mean(read.fun(which="V-theta2-unadj", n=sample.size)=="A3")
        }))
    
    ranked.theta.censoring.unadj <- data.frame(censoring="unadjusted weights (c)",
                                               sample.size=c(100, 200, 500, 1000, 1500, 2000),
                                               treatment=c(rep("A[1]~(with~direct~effect~on~event~of~interest)", length(ranked.theta1.A1)),
                                                           rep("A[2]~(with~effect~only~on~the~competing~risk~event)", length(ranked.theta1.A1)),
                                                           rep("A[3]~(with~no~effect)", length(ranked.theta1.A1))),
                                               theta=c(rep("bar(theta)[net]", 3*length(ranked.theta1.A1)),
                                                       rep("bar(theta)[crude]", 3*length(ranked.theta2.A1))),
                                               rank=c(ranked.theta1.A1,
                                                      ranked.theta1.A2,
                                                      ranked.theta1.A3,
                                                      ranked.theta2.A1,
                                                      ranked.theta2.A2,
                                                      ranked.theta2.A3))

    
    ranked.theta <- rbind(ranked.theta, ranked.theta.censoring.max, ranked.theta.censoring.unadj)

    (p.ranking <- ggplot(ranked.theta[ranked.theta$censoring!=
                                      "adjusted weights (a)",]) +
         #geom_point(aes(x=sample.size, y=rank, shape=theta, col=treatment)) +
         #geom_line(aes(x=sample.size, y=rank, lty=theta, col=treatment)) +
         geom_point(aes(x=sample.size, y=rank, shape=theta, col=treatment), size=2) +
         geom_line(aes(x=sample.size, y=rank, lty=treatment, col=theta)) +
         theme_bw() + ylab("fraction of times ranked no. 1") + xlab("sample size (n)") + 
         scale_fill_manual("", values=1:3, labels=parse_format()) +
         scale_linetype_manual("", values=1:3, labels=parse_format()) +
         scale_shape_manual("", values=1:3, labels=parse_format()) +
         #scale_color_manual("", values=hue_pal()(3)[1:3], labels=parse_format()) +
         scale_color_manual("", values=rep("black", 5), labels=parse_format()) +
         theme(axis.text=element_text(size=7),
               axis.title=element_text(size=8),
               title=element_text(size=8),
               legend.text=element_text(size=8),
               strip.text.x=element_text(size=8),
               strip.background=element_rect(fill="white", colour="black"),
               legend.position="bottom") +
         ggtitle("Ranking effectiveness (fraction of times ranked no. 1)") +
         #facet_wrap(. ~ censoring) +#, labeller=label_parsed) +
         facet_wrap(. ~ censoring)+#, labeller=label_parsed) +
         #guides(colour = guide_legend(nrow=5)))
         guides(colour = FALSE) +
         guides(lty = guide_legend(nrow=5)))

    ggsave(paste0("~/research/phd/random-forest/worg/figures/",
                  "plot-ranking-sample-size", ".pdf"),
           p.ranking,
           width=6.5, height=5)

    ggsave(paste0("~/Dropbox/PhD/org-helene/figures/",
                  "plot-ranking-sample-size", ".pdf"),
           p.ranking,
           width=6.5, height=5) 
    
    (p.ranking.b <- ggplot(ranked.theta[ranked.theta$censoring==
                                        "adjusted weights (b)",]) +
         #geom_point(aes(x=sample.size, y=rank, shape=theta, col=treatment)) +
         #geom_line(aes(x=sample.size, y=rank, lty=theta, col=treatment)) +
         geom_point(aes(x=sample.size, y=rank, shape=theta, col=treatment), size=2) +
         geom_line(aes(x=sample.size, y=rank, lty=treatment, col=theta)) +
         theme_bw() + ylab("fraction of times ranked no. 1") + xlab("sample size (n)") + 
         scale_fill_manual("", values=1:3, labels=parse_format()) +
         scale_linetype_manual("", values=1:3, labels=parse_format()) +
         scale_shape_manual("", values=1:3, labels=parse_format()) +
         #scale_color_manual("", values=hue_pal()(3)[1:3], labels=parse_format()) +
         scale_color_manual("", values=rep("black", 5), labels=parse_format()) +
         theme(axis.text=element_text(size=7),
               axis.title=element_text(size=8),
               title=element_text(size=8),
               legend.text=element_text(size=8),
               strip.text.x=element_text(size=8),
               strip.background=element_rect(fill="white", colour="black"),
               legend.position="right") +
         ggtitle("Ranking effectiveness (fraction of times ranked no. 1)") +
         #facet_wrap(. ~ censoring) +#, labeller=label_parsed) +
         #facet_wrap(. ~ censoring)+#, labeller=label_parsed) +
         #guides(colour = guide_legend(nrow=5)))
         guides(colour = FALSE) +
         ylim(c(0,1)) + 
         guides(lty = guide_legend(nrow=5)))

    ggsave(paste0("~/research/phd/random-forest/worg/figures/",
                  "plot-ranking-sample-size-onlyb", ".pdf"),
           p.ranking.b,
           width=8.5*0.9, height=4*0.8)

    (p.ranking.c <- ggplot(ranked.theta[ranked.theta$censoring==
                                        "unadjusted weights (c)",]) +
         #geom_point(aes(x=sample.size, y=rank, shape=theta, col=treatment)) +
         #geom_line(aes(x=sample.size, y=rank, lty=theta, col=treatment)) +
         geom_point(aes(x=sample.size, y=rank, shape=theta, col=treatment), size=2) +
         geom_line(aes(x=sample.size, y=rank, lty=treatment, col=theta)) +
         theme_bw() + ylab("fraction of times ranked no. 1") + xlab("sample size (n)") + 
         scale_fill_manual("", values=1:3, labels=parse_format()) +
         scale_linetype_manual("", values=1:3, labels=parse_format()) +
         scale_shape_manual("", values=1:3, labels=parse_format()) +
         #scale_color_manual("", values=hue_pal()(3)[1:3], labels=parse_format()) +
         scale_color_manual("", values=rep("black", 5), labels=parse_format()) +
         theme(axis.text=element_text(size=7),
               axis.title=element_text(size=8),
               title=element_text(size=8),
               legend.text=element_text(size=8),
               strip.text.x=element_text(size=8),
               strip.background=element_rect(fill="white", colour="black"),
               legend.position="right") +
         ggtitle("Ranking effectiveness (fraction of times ranked no. 1)") +
         #facet_wrap(. ~ censoring) +#, labeller=label_parsed) +
         #facet_wrap(. ~ censoring)+#, labeller=label_parsed) +
         #guides(colour = guide_legend(nrow=5)))
         ylim(c(0,1)) + 
         guides(colour = FALSE) +
         guides(lty = guide_legend(nrow=5)))

    ggsave(paste0("~/research/phd/random-forest/worg/figures/",
                  "plot-ranking-sample-size-onlyc", ".pdf"),
           p.ranking.c,
           width=8.5*0.9, height=4*0.8)
    
  
}



if (FALSE) { ### coverage as a function of sample size

    read.fun <- function(which="V-theta1-adj", n=1000, M=500, which.A=1,
                         new="-censoring",
                         censoring="medium") {
        ATE.list <- readRDS(file=paste0("output/ATE-list-fast",
                                        "-", which, "-fast-n",n,"-M", M, new, censoring, ".rds"))
        return(data.frame(which=paste0("A[", which.A, "]"), sample.size=round(n), coverage=round(coverage.fun(ATE.list)[which.A], 3)))
    }

    (coverage.theta1.A1 <-
         do.call("rbind", lapply(c(#100,
                              100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                  read.fun(n=sample.size, which.A=1)
                              })))

    (coverage.theta1.A2 <-
         do.call("rbind", lapply(c(#100,
                              100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                  read.fun(n=sample.size, which.A=2)
                              })))

    (coverage.theta1.A3 <-
         do.call("rbind", lapply(c(#100,
                              100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                  read.fun(n=sample.size, which.A=3)
                              })))

    (coverage.theta1.A1 <-
         do.call("rbind", lapply(c(#100,
                              100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                  read.fun(n=sample.size, which.A=1,
                                           censoring="max")
                              })))

    (coverage.theta1.A2 <-
         do.call("rbind", lapply(c(#100,
                              100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                  read.fun(n=sample.size, which.A=2,
                                           censoring="max")
                              })))

    (coverage.theta1.A3 <-
         do.call("rbind", lapply(c(#100,
                              100, 200, 500, 1000, 1500, 2000), function(sample.size) {
                                  read.fun(n=sample.size, which.A=3,
                                           censoring="max")
                              })))


    ## coverage.out <- lapply(c("I-theta2", "II-theta2-unadj", "II-theta2-adj", "III-theta1", "III-theta2", "IV-theta1-unadj", "IV-theta1-adj", "IV-theta2"
    ##                          ), function(which) {
    ##                              do.call("rbind", lapply(c(100, 100, 200, 500, 1000, 1500), function(sample.size) read.fun(which=which, n=sample.size)))
    ##                          })

}



if (FALSE) {
    
    ranked.1.n1000 <- unlist(lapply(c("I-theta2", "II-theta2-unadj", "II-theta2-adj", "III-theta1", "III-theta2", "IV-theta1-unadj", "IV-theta1-adj", "IV-theta2"
                                      ), function(which) {
                                          mean(read.fun(which=which)=="A1")
                                      }
                                    ))

    ranked.1.n100 <- unlist(lapply(c("I-theta2", "II-theta2-unadj", "II-theta2-adj", "III-theta1", "III-theta2", "IV-theta1-unadj", "IV-theta1-adj", "IV-theta2"
                                     ), function(which) {
                                         mean(read.fun(n=100, which=which)=="A1")
                                     }
                                   ))

    ranked.1.n200 <- unlist(lapply(c("I-theta2", "II-theta2-unadj", "II-theta2-adj", "III-theta1", "III-theta2", "IV-theta1-unadj", "IV-theta1-adj", "IV-theta2"
                                     ), function(which) {
                                         mean(read.fun(n=200, which=which)=="A1")
                                     }
                                   ))
    
    ranked.1.n500 <- unlist(lapply(c("I-theta2", "II-theta2-unadj", "II-theta2-adj", "III-theta1", "III-theta2", "IV-theta1-unadj", "IV-theta1-adj", "IV-theta2"
                                     ), function(which) {
                                         mean(read.fun(n=500, which=which)=="A1")
                                     }
                                   ))

    ranked.1.n1500 <- unlist(lapply(c("I-theta2", "II-theta2-unadj", "II-theta2-adj", "III-theta1", "III-theta2", "IV-theta1-unadj", "IV-theta1-adj", "IV-theta2"
                                     ), function(which) {
                                         mean(read.fun(n=1000, which=which)=="A1")
                                     }
                                   ))

}





if (FALSE) { ### coverage as a function of sample size

    read.fun <- function(which="I-theta2", n=1000, M=1000, which.A=1) {
        ATE.list <- readRDS(file=paste0("output/ATE-list-fast",
                                        "-", which, "-fast-n",n,"-M",M,new, censoring, ".rds"))
        return(data.frame(which=which, sample.size=n, coverage=round(coverage.fun(ATE.list)[which.A], 3)))
    }

    coverage.out <- lapply(c("I-theta2", "II-theta2-unadj", "II-theta2-adj", "III-theta1", "III-theta2", "IV-theta1-unadj", "IV-theta1-adj", "IV-theta2"
                             ), function(which) {
                                 do.call("rbind", lapply(c(100, 100, 200, 500, 1000, 1500), function(sample.size) read.fun(which=which, n=sample.size)))
                             })

    tab1 <- do.call("rbind", lapply(coverage.out, function(cov) {
        out <- rbind(t(cov)[3,])
        rownames(out) <- gsub("theta1", "$\\bar{\\theta}_1$",
                              gsub("theta2", "$\\bar{\\theta}_2$",
                                   gsub("V-", "V: ",
                                        gsub("I-", "I: ",
                                             gsub("-adj", " (adjusted weights)",
                                                  gsub("-unadj", " (unadjusted weights)", as.character(cov[1,1])))))))
        colnames(out) <- paste0("$n=", cov[, 2], "$")
        return(out)
    }))

    library(xtable)

    x <- xtable(tab1,
                align=paste0(rep("l", ncol(tab1)+1), collapse=""),
                caption=paste0("Coverage for estimation of the effect of $A_1$ as a function of sample size"),
                label=paste0("table:overview:coverage"))

    print(x, hline.after=c(0,nrow(x)),
          include.colnames=TRUE,
          include.rownames=TRUE,
          booktabs=TRUE,
          sanitize.text.function=function(x){x})  

}
