### misspecified_targets.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May  5 2022 (11:10) 
## Version: 
## Last-Updated: Jun  3 2022 (07:32) 
##           By: Thomas Alexander Gerds
##     Update #: 11
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
# misspecified censoring
misspecified_fixed <- list(horizon = 5,
                          event.times = c("T1","T2"),
                          treatments = c("A1" = .4,"A2" = .3, "A3" = .3,"A4" = .4,"A5" = .5,"A6" = .2,"A7" = .7,"A8" = .8,"A9" = .9,"A10" = .1),
                          binary.covariates = c("X1" = .1,"X2" = .2,"X3" = .3,"X4" = .4,"X5" = .5),
                          normal.covariates = paste0("X",6:7),
                          quadratic.covariates = "X6",
                          formula.list =  list(T1 ~ f(A4,.3) + f(A5,.7) + f(X1,1) + f(X2,.3) + f(X6,.3)+f(X6_2,.4),
                                               T2  ~ f(A5, -.3)+ f(X1,-.1) + f(X2,.6) + f(X6,.1),
                                               C ~ f(A1,0),
                                               A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2)))
misspecified_fixed_target <- tar_target(MISSPECIFIED_FIXED,misspecified_fixed,deployment = "main")
misspecified_varying <- data.table::CJ(A1_T1 = c(.7), #c(1,1.25,1.5)
                                      A1_T2 = 1,
                                      A2_T1 = 1,
                                      A2_T2 = 1.25,
                                      scale.censored = 1/25,# c(-Inf,1/40,1/25),
                                      sample.size = c(5000))
misspecified_varying_target <- tar_target(MISSPECIFIED_VARYING,misspecified_varying,deployment = "main")
misspecified_estimators <- rbindlist(list(data.table::CJ(num.trees = c(100),
                                                        net = c(FALSE),
                                                        method = c("causal_forest"),
                                                        weighter = c("ranger")),
                                         data.table::CJ(num.trees = c(100),
                                                        net = c(FALSE),
                                                        method = c("naive","CSC","FGR"),
                                                        weighter = c("km"))))
misspecified_estimators_target <- tar_target(MISSPECIFIED_ESTIMATORS,misspecified_estimators,deployment = "main")

# ---------------------------------------------------------------------
# Calculation of true ATE values
# ---------------------------------------------------------------------
#   crude: real world with competing risks
#     net: hypothetical world where competing risks are eliminated
misspecified_truth_varying <- tar_map(
    # truth is affected by % censored and true effect
    values = misspecified_varying[,.(scale.censored,A1_T1,A1_T2,A2_T1,A2_T2)],
    tar_target(MISSPECIFIED_VARYING_TRUTH,{
        lavaModel;
        simulateData;
        y = theTruth(setting = MISSPECIFIED_FIXED,
                     A1_T1 = A1_T1,
                     A1_T2 = A1_T2,
                     A2_T1 = A2_T1,
                     A2_T2 = A2_T2,
                     scale.censored = scale.censored,
                     horizon = MISSPECIFIED_FIXED$horizon,
                     n = 1000000)
        y = cbind(y, A1_T1 = A1_T1,A1_T2 = A1_T2, A2_T1 = A2_T1,A2_T2 = A2_T2)
        y},
        deployment = "worker"))

misspecified_truth <- tar_combine(MISSPECIFIED_TRUTH,
                                 misspecified_truth_varying)
# ---------------------------------------------------------------------
# Estimation of ATE 
# ---------------------------------------------------------------------
misspecified_estimates <- tar_map(
    # outer map generates data under varying parameters of the data generating model
    values = misspecified_varying[,.(A1_T1,A1_T2,A2_T1,A2_T2,scale.censored,sample.size)],
    unlist = FALSE,
    tar_target(MISSPECIFIED_SIM_DATA,{
        lavaModel
        simulateData(setting = MISSPECIFIED_FIXED,
                     A1_T1 = A1_T1,
                     A1_T2 = A1_T2,
                     A2_T1 = A2_T1,
                     A2_T2 = A2_T2,
                     n = sample.size,
                     scale.censored = scale.censored,
                     keep.latent = FALSE)},
        deployment = "worker",
        pattern = map(REPETITIONS)),
    # inner map estimates ATE under varying hyper-parameters of the estimator
    inner_misspecified <- tar_map(
        values = misspecified_estimators[,.(net,num.trees,method,weighter)],
        unlist = FALSE,
        tar_target(MISSPECIFIED_ESTIMATE,{
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
                              data=MISSPECIFIED_SIM_DATA,
                              times=MISSPECIFIED_FIXED$horizon,
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
        pattern = map(MISSPECIFIED_SIM_DATA))
    )
)
# combine
misspecified_ate <- tar_combine(MISSPECIFIED_ESTIMATE_ATE,{
    # the first element are the data
    misspecified_estimates[-1]
})
# ---------------------------------------------------------------------
# Summarize performance of estimators against true parameter values
# ---------------------------------------------------------------------
misspecified_results <-
    tar_target(MISSPECIFIED_RESULTS,
               summarizePerformance(truth = MISSPECIFIED_TRUTH,
                                    estimate = MISSPECIFIED_ESTIMATE_ATE),
               deployment = "main")

misspecified_boxplots <- tar_target(MISSPECIFIED_BOXPLOTS,{
    e <- MISSPECIFIED_ESTIMATE_ATE[intervene == "A1"]
    t <- MISSPECIFIED_TRUTH
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
    g+facet_grid(~A1_T1)
})
misspecified_coverage <- tar_target(MISSPECIFIED_COVERAGE,{
    x <- MISSPECIFIED_RESULTS
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
misspecified_plot_setting <- tar_map(
    values = misspecified_varying[,.(scale.censored,A1_T1,A1_T2,A2_T1,A2_T2)],
    tar_target(PLOT_MISSPECIFIED,
               plotTruth(setting = MISSPECIFIED_FIXED,
                         A1_T1 = A1_T1,
                         A1_T2 = A1_T2,
                         A2_T1 = A2_T1,
                         A2_T2 = A2_T2,
                         scale.censored = scale.censored,
                         horizon = MISSPECIFIED_FIXED$horizon,
                         n = 10000)))




######################################################################
### misspecified_targets.R ends here
