try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
try(setwd("/maps/projects/biostat01/people/grb615/research/SoftWare/grfCausalSearch"),silent=TRUE)
library(targets)
## e <- tar_read(ESTIMATE_ATE)
## e[,n:=factor(n,levels=c("500","1000","2000","5000"),labels=c("500","1000","2000","5000"))]
## e[,net:=factor(net,levels=c("0","1"),labels=c("0","1"))]
## e <- e[theme == "ranking"]
## e2 <- e[intervene == "A2"&A2_T2 == 2]
## ggplot(e2,aes(y = ate,group = net,color = net))+geom_boxplot()+facet_grid(~n)
library(parallel)
library(targets)
library(tarchetypes)
library(ranger)
library(survival)
source("./setting/simulation_targets.R")
for(f in list.files("R",".R$",full.names=TRUE)){source(f)}
for(f in list.files("functions",".R$",full.names=TRUE)){source(f)}
this <- varying_ranking[A2_T2 == 2&sample.size == 1000&net == TRUE]
this$A2_T2 <- .1
this$horizon <- 5
library(data.table)

set.seed(8159)
REPETITIONS <- 1:20
MC <- 5
## formula_ranking <- list(T1 ~ f(X1,0) + f(X2,0) + f(X6,0),T2  ~ f(X1,0) + f(X2,0) + f(X6,0),C ~ f(A1,0),A1 ~ f(X1,0) + f(X6,0) + f(A7,0))
formula_ranking <- list(T1 ~ f(X1,1) + f(X2,.3) + f(X6,.3),T2  ~ f(X1,-.1) + f(X2,.6) + f(X6,.1),C ~ f(A1,0),A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2))
## formula_ranking <- list(T1 ~ f(X1,1) + f(X2,.3) + f(X6,.3),T2  ~ f(X1,-.1) + f(X2,.3) + f(X6,-0.1),C ~ f(A1,0),A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2))
this$scale.censored <- -Inf
this$num.trees <- 50
Y <- do.call("rbind",lapply(c(FALSE),function(net){
    do.call("rbind",mclapply(REPETITIONS,function(b){
        print(b)
        set = fixed
        set$formula.list = formula_ranking
        tt = theTruth(setting = set,A1_T1 = this$A1_T1,A1_T2 = this$A1_T2,A2_T1 = this$A2_T1,A2_T2 = this$A2_T2,horizon = this$horizon,scale.censored = this$scale.censored,B = 1,cores = 1,n = 100000)[intervene == "A2"&cause == 1]
        simulated_data <- simulateData(setting = set,A1_T1 = this$A1_T1,A1_T2 = this$A1_T2,A2_T1 = this$A2_T1,A2_T2 = this$A2_T2,n = this$sample.size,scale.censored = this$scale.censored,keep.latent = FALSE)
        ## ff = Hist(time,event)~intervene(A1)+intervene(A2)+intervene(A3)+intervene(A4)+intervene(A5)+intervene(A6)+intervene(A7)+intervene(A8)+intervene(A9)+intervene(A10)+X1+X2+X3+X4+X5+X6+X7
        ff = Hist(time,event)~A1+intervene(A2)+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7
        ## w1 = weighter(formula = Hist(time,event)~A1+A2,data = simulated_data,CR.as.censoring = 1,times = 5,method = "km")
        ## w2 = weighter(formula = Hist(time,event)~A1+A2,data = simulated_data,CR.as.censoring = 1,times = 5,method = "ranger")
        ## cbind(w1[w1 != 0],w2[w2 != 0])
        ## plot(w1[w1 != 0],w2[w2 != 0])
        ## cbind(w1,w2)
        x1 <- causalhunter(formula=ff,method = "causal_forest",weighter="km",args.weight = list(num.trees = this$num.trees),num.trees=this$num.trees,CR.as.censoring = net,data=simulated_data,times=this$horizon,formula.weight = Hist(time,event)~A1+A2)
        x2 <- causalhunter(formula=ff,method = "causal_forest",weighter="ranger",args.weight = list(num.trees = this$num.trees,alpha = 0.05,mtry = 17),fit.separate = TRUE,num.trees=this$num.trees,CR.as.censoring = net,data=simulated_data,times=this$horizon,formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)
        x3 <- causalhunter(formula=ff,method = "causal_forest",weighter="ranger",args.weight = list(num.trees = 50,alpha = 0.05,mtry = 17),fit.separate = TRUE,num.trees=this$num.trees,CR.as.censoring = net,data=simulated_data,times=this$horizon,formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)
        ## formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)
        ## x[,rank := rank(-abs(ate))]
        ## x <- cbind(x,data.table(n = this$sample.size,net = as.numeric(net),scale.censored = this$scale.censored,num.trees = this$num.trees,method = "causal_forest",A1_T1 = this$A1_T1,A1_T2 = this$A1_T2,A2_T1 = this$A2_T1,A2_T2 = this$A2_T2,formula = "this",theme = "ranking"))
        x <- cbind(x1[,data.table::data.table(intervene,km = ate,net = net)],x2[,data.table::data.table(ranger = ate)],x3[,data.table::data.table(ranger2 = ate)])
        ## if (net == 1) x = cbind(x,truth = tt[net == 1]$ate) else x = cbind(x,truth = tt[net == 0]$ate)
        x
    },mc.cores = MC))}))
y <- Y[,data.table::data.table(ranger,ranger2,km,net)]
y[,data.table::data.table(km = mean(km),ranger = mean(ranger),ranger2 = mean(ranger2)),by = net]



y1 <- Y[intervene == "A1",data.table::data.table(ate,net,n)]
y1[,mean(ate),by = net]

