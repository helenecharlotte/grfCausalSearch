library(grf)
library(data.table)
## library(rlist)
library(survival)
library(riskRegression)
library(ggplot2)
library(gridExtra)

source("sim-data.R")
source("hunt-fun.R")
source("hunt-fun-fast.R")


if (FALSE) {
    library("devtools")
    setwd("/home/helene/research/phd/random-forest/grf-causal-search-master")
    install_github("helenecharlotte/grfCausalSearch")
    devtools::document()
    devtools::build()
    devtools::install()
    devtools::check()
}


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
