#-------------------------------------------------------------------------------------------#
## set working directory
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){
    setwd("/home/ifsv/jhl781/research/phd/random-forest/r-code/causal-search-dst/")
} else {
    setwd("~/research/phd/random-forest/r-code/causal-search-dst/")
}


library(grf)
library(data.table)
library(survival)
library(riskRegression)

source("sim-data.R")
source("hunt-fun-fast.R")


#--------------------- test with simulations ------------------------#'


n <- 1000 # number of subjects
M <- 10 # number of simulations
n.A <- 10 # number of treatments

#-- compute true values for all treatments and both theta1 and theta2:

psi0.list.I.2 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=2, CR=1, which.A=a)}))
psi0.list.II.2 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=2, CR=1, which.A=a, cens.A=10)}))
psi0.list.III.1 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=1, CR=2, which.A=a,
             form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                      as.numeric(X[, 3])*0.33 +
                                      0*A[, 2])}))
psi0.list.III.2 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=2, CR=2, which.A=a,
             form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                      as.numeric(X[, 3])*0.33 +
                                      0*A[, 2])}))
psi0.list.IV.1 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=1, CR=2, which.A=a)}))
psi0.list.IV.2 <- unlist(lapply(1:n.A, function(a) {
    sim.data(1e6, compute.psi=2, CR=2, which.A=a)}))


names(psi0.list.I.2) <- names(psi0.list.II.2) <- paste0("A", 1:10)
names(psi0.list.III.1) <- names(psi0.list.III.2) <- paste0("A", 1:10)
names(psi0.list.IV.1) <- names(psi0.list.IV.2) <- paste0("A", 1:10)

#-- repeat simulations:

ATE.list.I.2 <- list(true.values=psi0.list.I.2)
ATE.list.II.2.unadj <- list(true.values=psi0.list.II.2)
ATE.list.II.2.adj <- list(true.values=psi0.list.II.2)

ATE.list.III.1 <- list(true.values=psi0.list.III.1)
ATE.list.III.2 <- list(true.values=psi0.list.III.2)

ATE.list.IV.1.unadj <- list(true.values=psi0.list.IV.1)
ATE.list.IV.1.adj <- list(true.values=psi0.list.IV.1)
ATE.list.IV.2 <- list(true.values=psi0.list.IV.2)

for (m in 1:M) {

    #---- scenario I:
    dt <- sim.data(n, CR=1, seed=201010120+m)

    # estimate theta2:
    ATE.list.I.2[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                         X.vars=paste0("X", 1:10))

    #---- scenario II:
    dt <- sim.data(n, CR=1, cens.A=10, seed=201010120+m)

    # estimate theta2:
    ATE.list.II.2.unadj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                                X.vars=paste0("X", 1:10))
    ATE.list.II.2.adj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                              X.vars=paste0("X", 1:10),
                                              stratify.CR="A10")

    #---- scenario III:
    dt <- sim.data(n, CR=2, form.T2 = function(X, A)
        0.1 - as.numeric(X[, 2])*0.4 -
        as.numeric(X[, 3])*0.33 +
        0*A[, 2],
        seed=201010120+m)

    # estimate theta1:
    ATE.list.III.1[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                           X.vars=paste0("X", 1:10))

    # estimate theta2:
    ATE.list.III.2[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                           X.vars=paste0("X", 1:10))

    #---- scenario IV:
    dt <- sim.data(n, CR=2, seed=201010120+m)

    # estimate theta1:
    ATE.list.IV.1.unadj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                                X.vars=paste0("X", 1:10))
    ATE.list.IV.1.adj[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=TRUE,
                                              X.vars=paste0("X", 1:10),
                                              stratify.CR="A2")
    
    # estimate theta2:
    ATE.list.IV.2[[m+1]] <- hunt.fun.fast(dt, CR.as.censoring=FALSE,
                                          X.vars=paste0("X", 1:10))

}

