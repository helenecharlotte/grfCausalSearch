#-------------------------------------------------------------------------------------------#
## set working directory
#-------------------------------------------------------------------------------------------#

library(grf)
library(data.table)
library(survival)
library(riskRegression)

source("R/sim-data.R")
source("R/hunt-fun-fast.R")


#--------------------- test with simulations ------------------------#'


n <- 1000 # number of subjects
M <- 10 # number of simulations
n.A <- 10 # number of treatments

censoring <- "medium" #"medium" # "max", "min"

if (censoring=="medium") {
    C.shape <- 0.4
} else if (censoring=="max") {
    C.shape <- 0.2
} else if (censoring=="min") {
    C.shape <- 0.6
}


#-- compute true values for all treatments and both theta1 and theta2:

psi0.list.1 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=1, CR=2, which.A=a,
             form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                      as.numeric(X[, 1])*0.33 +
                                      1.5*A[, 2])}))
psi0.list.2 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=2, CR=2, which.A=a,
             form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                      as.numeric(X[, 1])*0.33 +
                                      1.5*A[, 2])}))

names(psi0.list.1) <- names(psi0.list.2) <- paste0("A", 1:10)

#-- repeat simulations:

ATE.list.1.unadj <- list(true.values=psi0.list.1)
ATE.list.1.adjA <- list(true.values=psi0.list.1)
ATE.list.1.adj <- list(true.values=psi0.list.1)
ATE.list.2.unadj <- list(true.values=psi0.list.2)
ATE.list.2.adjA <- list(true.values=psi0.list.2)
ATE.list.2.adj <- list(true.values=psi0.list.2)

nevents <- list()

for (m in 1:M) {

    dt <- sim.data(n, CR=2, seed=201010120+m,
                   C.shape=C.shape,
                   form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                            as.numeric(X[, 1])*0.33 +
                                            1.5*A[, 2])

    nevents[[m]] <- table(dt[time<=0.5, delta])

    # estimate theta1:
    ATE.list.1.unadj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                             X.vars=paste0("X", 1:10))
    ATE.list.1.adjA[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                            X.vars=paste0("X", 1:10),
                                            stratify.CR="A2")
    ATE.list.1.adj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                           X.vars=paste0("X", 1:10),
                                           stratify.CR="A2+X1+X2")

    # estimate theta2:
    ATE.list.2.unadj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                               X.vars=paste0("X", 1:10))
    ATE.list.2.adjA[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                              X.vars=paste0("X", 1:10),
                                              stratify.CR="A2")
    ATE.list.2.adj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                             X.vars=paste0("X", 1:10),
                                             stratify.CR="A2+X1+X2")

}

