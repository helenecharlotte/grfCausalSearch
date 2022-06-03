### net_targets.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May  5 2022 (11:08) 
## Version: 
## Last-Updated: May 23 2022 (08:43) 
##           By: Thomas Alexander Gerds
##     Update #: 23
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# net censoring
library(data.table)
net_fixed <- list(horizon = 5,
                  event.times = c("T1","T2"),
                  treatments = c("A1" = .4,"A2" = .3, "A3" = .3,"A4" = .4,"A5" = .5,"A6" = .2,"A7" = .7,"A8" = .8,"A9" = .9,"A10" = .1),
                  binary.covariates = c("X1" = .1,"X2" = .2,"X3" = .3,"X4" = .4,"X5" = .5),
                  normal.covariates = paste0("X",6:7),
                  formula.list =  list(T1 ~ f(A4,.3) + f(A5,.7) + f(X1,1) + f(X2,.3) + f(X6,-.5),
                                       T2  ~ f(A5, -.3)+ f(X1,-.1) + f(X2,.6) + f(X6,.1),
                                       C ~ f(A1,0),
                                       A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2)))
net_fixed_target <- tar_target(NET_FIXED,net_fixed,deployment = "main")
net_varying <- data.table::CJ(A1_T1 = c(.8,1,1.25,1.5),
                              A1_T2 = c(.8,1,1.25),
                              A2_T1 = 1,
                              A2_T2 = 1.25,
                              scale.censored = 1/40,
                              sample.size = c(500,1000,2000))
net_varying_target <- tar_target(NET_VARYING,net_varying,deployment = "main")
net_estimators <- rbindlist(list(data.table::CJ(num.trees = c(50,100,150),
                                                net = TRUE,
                                                method = c("causal_forest"),
                                                weighter = c("ranger"))))
net_estimators_target <- tar_target(NET_ESTIMATORS,net_estimators,deployment = "main")

# ---------------------------------------------------------------------
# Calculation of true ATE values
# ---------------------------------------------------------------------
#   crude: real world with competing risks
#     net: hypothetical world where competing risks are eliminated
net_truth_varying <- tar_map(
    # truth is affected by % censored and true effect
    values = net_varying[,.(scale.censored,A1_T1,A1_T2,A2_T1,A2_T2)],
    tar_target(NET_VARYING_TRUTH,{
        lavaModel;
        simulateData;
        y = theTruth(setting = NET_FIXED,
                     A1_T1 = A1_T1,
                     A1_T2 = A1_T2,
                     A2_T1 = A2_T1,
                     A2_T2 = A2_T2,
                     scale.censored = scale.censored,
                     horizon = NET_FIXED$horizon,
                     n = 1000000)
        y = cbind(y, A1_T1 = A1_T1,A1_T2 = A1_T2, A2_T1 = A2_T1,A2_T2 = A2_T2)
        y},
        deployment = "worker"))

net_truth <- tar_combine(NET_TRUTH,net_truth_varying)
# ---------------------------------------------------------------------
# Estimation of ATE 
# ---------------------------------------------------------------------
net_estimates <- tar_map(
    # outer map generates data under varying parameters of the data generating model
    values = net_varying[,.(A1_T1,A1_T2,A2_T1,A2_T2,scale.censored,sample.size)],
    unlist = FALSE,
    tar_target(NET_SIM_DATA,{
        lavaModel
        simulateData(setting = NET_FIXED,
                     A1_T1 = A1_T1,
                     A1_T2 = A1_T2,
                     A2_T1 = A2_T1,
                     A2_T2 = A2_T2,
                     n = sample.size,
                     scale.censored = scale.censored,
                     keep.latent = FALSE)},
        deployment = "main",
        pattern = map(REPETITIONS)),
    # inner map estimates ATE under varying hyper-parameters of the estimator
    inner_net <- tar_map(
        values = net_estimators[,.(net,num.trees,method,weighter)],
        unlist = FALSE,
        tar_target(ESTIMATE,{
            if(weighter == "km")
                ff <- Hist(time,event)~A1+A2
            else
                ff <- Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7
            x <- causalhunter(formula=Hist(time,event)~intervene(A1)+intervene(A2)+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7,
                              method = method,
                              weighter=weighter,
                              args.weight = list(num.trees = num.trees),
                              num.trees=num.trees,
                              CR.as.censoring = net,
                              data=NET_SIM_DATA,
                              times=NET_FIXED$horizon,
                              formula.weight = ff)
            x <- cbind(x,data.table(n = sample.size,
                                    net = as.numeric(net),
                                    scale.censored = scale.censored,
                                    num.trees = num.trees,
                                    method = method,
                                    weighter = weighter,
                                    A1_T1 = A1_T1,
                                    A1_T2 = A1_T2,
                                    A2_T1 = A2_T1,
                                    A2_T2 = A2_T2))
            x
        },
        deployment = "worker",
        pattern = map(NET_SIM_DATA))
    )
)
# combine
net_ate <- tar_combine(NET_ESTIMATE_ATE,{
    # the first element are the data
    net_estimates[-1]
})
# ---------------------------------------------------------------------
# Summarize performance of estimators against true parameter values
# ---------------------------------------------------------------------
net_results <-
    tar_target(NET_RESULTS,
               summarizePerformance(truth = NET_TRUTH,
                                    estimate = NET_ESTIMATE_ATE),
               deployment = "main")

net_boxplots <- tar_target(NET_BOXPLOTS,{
    e <- NET_ESTIMATE_ATE[intervene == "A1"]
    t <- NET_TRUTH
    setkeyv(e,c("net","intervene","A1_T1","A1_T2","A2_T1","A2_T2"))
    setkeyv(t,c("net","intervene","A1_T1","A1_T2","A2_T1","A2_T2"))
    t.ate = t[cause == 1,.(true.ate = mean(ate)),keyby = c("net","intervene","A1_T1","A1_T2","A2_T1","A2_T2")]
    e = merge(e,t.ate,by = c("net","intervene","A1_T1","A1_T2","A2_T1","A2_T2"),all.y = FALSE)
    e[,n:=factor(n)]
    e[,net:=factor(net)]
    e[,num.trees:=factor(num.trees)]
    e[,A1_T1:=factor(A1_T1)]
    e[,A1_T2:=factor(A1_T2)]
    e[,A2_T1:=factor(A2_T1)]
    e[,A2_T2:=factor(A2_T2)]
    e[method == "causal_forest",method := paste0(method,"_",weighter)]
    g <- ggplot(e,aes(x = method,y = ate))
    g <- g+geom_boxplot()
    g <- g+geom_hline(aes(yintercept = true.ate,color = "red"),data = e)
    g+facet_grid(n~A1_T1)
})
net_coverage <- tar_target(NET_COVERAGE,{
    x <- NET_RESULTS
    x[,censored.tau:=factor(censored.tau)]
    x[,n:=factor(n)]
    x[,net:=factor(net)]
    x[,num.trees:=factor(num.trees)]
    x[,A1_T1:=factor(A1_T1)]
    x[,A1_T2:=factor(A1_T2)]
    x[,A2_T1:=factor(A2_T1)]
    x[,A2_T2:=factor(A2_T2)]
    g <- ggplot(x[intervene == "A1"],aes(x = n,y = coverage,group = method,color = method))
    g <- g+geom_line()+facet_grid(A1_T1~censored.tau)
    g
})

# ---------------------------------------------------------------------
# Plot of setting
# ---------------------------------------------------------------------
net_plot_setting <- tar_map(
    values = net_varying[,.(scale.censored,A1_T1,A1_T2,A2_T1,A2_T2)],
    tar_target(PLOT_NET,
               plotTruth(setting = NET_FIXED,
                         A1_T1 = A1_T1,
                         A1_T2 = A1_T2,
                         A2_T1 = A2_T1,
                         A2_T2 = A2_T2,
                         scale.censored = scale.censored,
                         horizon = NET_FIXED$horizon,
                         n = 10000)))

######################################################################
### net_targets.R ends here
