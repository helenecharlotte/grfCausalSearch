### test.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 19 2022 (07:48) 
## Version: 
## Last-Updated: Jul 15 2022 (09:02) 
##           By: Thomas Alexander Gerds
##     Update #: 123
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
try(setwd("~/research/SoftWare/grfCausalSearch/"),silent=TRUE)
try(setwd("~/Dropbox/PhD/grfCausalSearch/"),silent=TRUE)
library(targets)
library(data.table)
library(splines)
library(rms)
library(Greg)
library(survival)
library(ranger)
library(survival)
library(parallel)
for(f in list.files("R",".R$",full.names=TRUE)){source(f)}
for(f in list.files("functions",".R$",full.names=TRUE)){source(f)}
fixed <- list(horizon = 5,
              event.times = c("T1","T2"),
              treatments = c("A1" = .4,"A2" = .3, "A3" = .3,"A4" = .4,"A5" = .5,"A6" = .6,"A7" = .7,"A8" = .8,"A9" = .9,"A10" = .1),
              binary.covariates = c("X1" = .1,"X2" = .2,"X3" = .3,"X4" = .4,"X5" = .5),
              normal.covariates = paste0("X",6:7),
              quadratic.covariates = "X6",
              scale = 1/100)
fixed$formula.list = list(T1 ~ f(A4,0.9) + f(A5,-.9) + f(X1,1) + f(X2,.3) + f(X6,.3) + f(X6_2,.4),
                          T2  ~ f(A3, -1.3)+ f(A4,1.3),
                          C ~ f(A1,0)+f(A3,0) + f(X2,0)+f(X6,0),
                          A1 ~ f(X1,0) + f(X6,0) + f(A7,0))

e <- tar_read(ESTIMATE_ATE)
e[,n:=factor(n,levels=c("500","1000","2000","5000"),labels=c("500","1000","2000","5000"))]
e[,net:=factor(net,levels=c("0","1"),labels=c("0","1"))]
e <- e[theme == "ranking"]
e2 <- e[intervene == "A2"&A2_T2 == 2]
ggplot(e2,aes(y = ate,group = net,color = net))+geom_boxplot()+facet_grid(~n)


## formula.list = list(T1 ~ f(A4,.3) + f(A5,.7) + f(X1,1) + f(X2,.3) + f(X6,.3) + f(X6_2,.4),
## formula.list = list(T1 ~ f(A4,.3) + f(A5,.7) + f(X1,1) + f(X2,.3) + f(X6,-.5),
formula_ranking <- list(T1 ~ f(X1,1) + f(X2,.3) + f(X6,.3),
                        T2  ~ f(X1,-.1) + f(X2,.6),
                        C ~ f(A1,0),
                        A1 ~ f(X1,-1) + f(X6,.7) + f(A7,.2))
fixed$formula.list <- formula_ranking



A1_T1 <- 1.25
A1_T2 <- 0.7
A2_T1 <- 1
A2_T2 <- 1.25
scale.censored <- 1/40
## fixed$formula.list <- formula_ranking
simulated_data <- simulateData(setting = fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,n = 2000,scale.censored = scale.censored,keep.latent = FALSE)
tt <- theTruth(setting = fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,scale.censored = scale.censored,horizon = fixed$horizon,B = 5,n = 100000,cores = 5)
tt[cause == 1 & intervene == "A1",.(intervene,ate = round(100*ate,2),net)]

ff = Hist(time,event)~intervene(A1)+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7

ff = Hist(time,event)~intervene(A1)+intervene(A2)+intervene(A3)+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7
X <- do.call("rbind",mclapply(1:50,function(b){
    print(b)
    simulated_data <- simulateData(setting = fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,n = 500,scale.censored = scale.censored,keep.latent = FALSE)
    xcrude <- causalhunter(formula=ff,
                           method = "causal_forest",
                           weighter="ranger",
                           args.weight = list(num.trees = 50,alpha = 0.05,mtry = 17),
                           num.trees=50,
                           CR.as.censoring = FALSE,
                           data=simulated_data,
                           times=fixed$horizon,
                           formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)
    xcrude
},mc.cores = 5))
X[,.(est = mean(ate),truth = tt[intervene == "A1"&net == 0&cause == 1]$ate)]


## publish(coxph(Surv(time,event == 1)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10,data = simulated_data))
## publish(coxph(Surv(time,event == 2)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10,data = simulated_data))
ff = Hist(time,event)~intervene(A1)+intervene(A2)+intervene(A3)+intervene(A4)+intervene(A5)+intervene(A6)+intervene(A7)+intervene(A8)+intervene(A9)+intervene(A10)+X1+X2+X3+X4+X5+X6+X7
x1 <- causalhunter(formula=ff,method = "causal_forest",weighter="ranger",args.weight = list(num.trees = 100),num.trees=100,CR.as.censoring = FALSE,data=simulated_data,times=fixed$horizon,formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)
x1[,rank := rank(-abs(ate))]
x1
## for(f in list.files("R",".R$",full.names=TRUE)){source(f)}

xnet <- causalhunter(formula=ff,method = "causal_forest",weighter="ranger",args.weight = list(num.trees = 50),num.trees=50,CR.as.censoring = TRUE,data=simulated_data,times=fixed$horizon,formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)
xnet[,rank := rank(-abs(ate))]
xnet

cbind(x1[,.(intervene,crude = round(100*ate,1),crude.rank = rank)],xnet[,.(intervene,net = round(100*ate,1),net.rank = rank)])


FIXED <- tar_read(FIXED)
## setting <- "formula1"
setting <- "formula_ranking"
set <- FIXED
set$formula.list = formula_settings[[setting]]
scale.censored <- 1/40
A1_T1 <- 1.25
A1_T2 <- 1
A2_T1 <- 1
A2_T2 <- 1.25
theTruth(setting = set,
         A1_T1 = A1_T1,
         A1_T2 = A1_T2,
         A2_T1 = A2_T1,
         A2_T2 = A2_T2,
         scale.censored = scale.censored,
         horizon = 5,
         B = 10,
         n = 100000,cores = 5)[cause == 1&net == 1&intervene == "A1"]



simulated_data <- simulateData(setting = fixed,A1_T1 = 1,A1_T2 = 1,A2_T1 = 1,A2_T2 = 1,n = 500,scale.censored = 1/25,keep.latent = FALSE)

x1 <- causalhunter(formula=ff,
                  method = "causal_forest",
                  weighter="ranger",
                  args.weight = list(num.trees = 100),
                  num.trees=100,
                  CR.as.censoring = FALSE,
                  data=simulated_data,
                  times=5,
                  formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7)


## formula.list = list(T1 ~ f(A4,0) + f(X1,0) + f(X2,0) + f(X6,0),
## T2  ~ f(A5, 0)+ f(X1,0) + f(X2,0) + f(X6,0),
## C ~ f(A1,0.3)+ f(A3,0.3) + f(A6,-.2) + f(X2,0.5)+f(X7,-.2),
## A1 ~ f(X1,0) + f(X6,0)))
scale.censored <- -Inf
## scale.censored <- 1/40
scale.censored <- 1/25
A1_T1 <- 1.3
A1_T2 <- 1
A2_T1 <- 1
A2_T2 <- 0.5
d <- simulateData(fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,n = 1000,scale.censored = scale.censored,keep.latent = FALSE)
fit.cph <- cph(Surv(time,event == 1)~rcs(X6),data=d)
fit.coxph <- cph(Surv(time,event == 1)~bs(X6,3),data=d)
plotHR(fit.cph, term = "X6", plot.bty = "o", xlim = c(-2, 2), xlab = "Age",ylim = c(.1,5.5))
plotHR(fit.coxph, term = "X6", plot.bty = "o", xlim = c(-2, 2), xlab = "Age")

plotTruth(setting = fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,scale.censored = scale.censored,horizon = 5,n = 1000)

y <- theTruth(setting = fixed,
              A1_T1 = A1_T1,
              A1_T2 = A1_T2,
              A2_T1 = A2_T1,
              A2_T2 = A2_T2,
              scale.censored = scale.censored,
              horizon = 3,
              n = 1000000)
print(y)


# compare with unadjusted KM weights
set.seed(81)
u <- mclapply(1:48,function(i){
    d <- simulateData(fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,n = 500,scale.censored = scale.censored,keep.latent = FALSE)
    c1 <- causalhunter(formula=Hist(time,event)~intervene(A1)+intervene(A2)+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7,weighter="km",num.tree=100,cause = 1,CR.as.censoring = 0,data=d,times=fixed$horizon,formula.weight = Hist(time,event)~1)
    c2s <- causalhunter(formula=Hist(time,event)~intervene(A1)+intervene(A2)+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7,weighter="ranger",num.trees=100,cause = 1,CR.as.censoring = 0,data=d,times=fixed$horizon,formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7,args.weight = list(num.trees = 100),fit.separate = TRUE)
    c(KM = c1$ate[1],maxstat = c2s$ate[1])
},mc.cores = 6)
boxplot(t(simplify2array(u)))
abline(h = y[intervene == "A1"&cause == 1&net == 0,ate],col = 2,lwd = 3)


# effect of sample size
set.seed(81)
u <- do.call("rbind",lapply(c(500,1000,1500,5000),function(n){
    u <- mclapply(1:48,function(i){
        d <- simulateData(fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,n = n,scale.censored = scale.censored,keep.latent = FALSE)
        c2s <- causalhunter(formula=Hist(time,event)~intervene(A1)+intervene(A2)+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7,weighter="ranger",num.trees=100,cause = 1,CR.as.censoring = 0,data=d,times=fixed$horizon,formula.weight = Hist(time,event)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7,args.weight = list(num.trees = 10),fit.separate = TRUE)
        c(maxstat = c2s$ate[1])
    } ,mc.cores = 6)
    c(n = n,ate = unlist(u))
}))
U <- t(u[,-1])
U[U[,2]>1,2] <- 0
boxplot(U)
abline(h = y[intervene == "A1"&cause == 1&net == 0,ate],col = 2,lwd = 3)


m = 20
ntree = 10
d2 <- simulateData(fixed,A1_T1 = A1_T1,A1_T2 = A1_T2,A2_T1 = A2_T1,A2_T2 = A2_T2,n = m,scale.censored = scale.censored,keep.latent = FALSE)
F = ranger(Surv(time,event == 0)~A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+X1+X2+X3+X4+X5+X6+X7,data = d,num.tree = ntree,keep.inbag = TRUE,min.node.size = 100,replace = FALSE)
x <- predict(F,data = d2,type = "terminalNodes")$predictions
inbag = 1*(do.call("cbind",F$inbag.counts) != 0)
y = predict(F,data = d, type = "terminalNodes")$predictions
Y = y*inbag
who = lapply(1:m,function(i){
    lapply(1:ntree,function(tree){which(x[i,tree] == Y[,tree])})
})
Who = lapply(1:m,function(i){do.call("c",who[[i]])})
G1 = lapply(1:m,function(i){Reduce("+",lapply(who[[i]],function(tree){
    predict(prodlim(Hist(time,event)~1,data = d[tree],reverse = TRUE,conf.int = FALSE),times = sort(unique(d$time)),type = "risk")
}))/ntree})
G2 = lapply(1:m,function(i){prodlim(Hist(time,event)~1,data = d[Who[[i]]],reverse = TRUE,conf.int = FALSE)})
lapply(G2,plot,add = TRUE)
G3 = prodlim(Hist(time,event)~1,data = d, reverse = TRUE,conf.int = FALSE)
plot(G3,col = 3,lwd = 3)
lapply(G2,plot,add = TRUE,col = "gray55")    
plot(G3,col = 3,lwd = 3,add = TRUE)
b = predictRisk(F,newdata = d2,times = sort(unique(d$time)))
nix = apply(b,1,function(g)lines(sort(unique(d$time)),g,col = 4))


### coverage in ranking settings
estimate <- tar_read("ESTIMATE_ATE")
truth <- tar_read("TRUTH")
tt <- truth[intervene == "A1"&formula == "formula_ranking"&cause == 1&net == 0]
ee <- estimate[intervene == "A1"&formula == "formula_ranking"&net == 0&n == 5000]
e02 <- ee[A2_T2 == 0.2,]
e2 <- ee[A2_T2 == 2,]
e02[,true.ate := tt[A2_T2 == 0.2,mean(ate)]]
e2[,true.ate := tt[A2_T2 == 2,mean(ate)]]
e02[,.(coverage = mean(lower <= true.ate & true.ate <= upper)),by = scale.censored]
e2[,.(coverage = mean(lower <= true.ate & true.ate <= upper)),by = scale.censored]

### coverage in crude settings
estimate <- tar_read("ESTIMATE_ATE")
truth <- tar_read("TRUTH")
tt <- truth[intervene == "A1"&formula == "formula1"&cause == 1&net == 0&horizon == 5]
ee <- estimate[intervene == "A1"&formula == "formula1"&net == 0&n == 5000&theme == "crude_effect"]
setkey(ee,A1_T1,A1_T2)
setkey(tt,A1_T1,A1_T2)
TT <- tt[,.(true.ate = ate,A1_T1,A1_T2)]
EE <- ee[,.(ate,se,lower,upper,A1_T1,A1_T2)]
u <- TT[EE]
u[,.(coverage = mean(lower <= true.ate & true.ate <= upper)),keyby = .(A1_T1,A1_T2)]
u[,mean(ate),by = A1_T2]
u[,table(true.ate)]

### coverage in censoring settings
estimate <- tar_read("ESTIMATE_ATE")
truth <- tar_read("TRUTH")
tt <- truth[intervene == "A1"&cause == 1&net == 0&horizon == 3]
ee <- estimate[intervene == "A1"&net == 0&n == 5000&theme == "censoring"][,.(ate,scale.censored,se,lower,upper,A1_T1,A1_T2)][]
setkey(ee,A1_T1,A1_T2)
setkey(tt,A1_T1,A1_T2)
TT <- tt[,.(true.ate = ate,A1_T1,A1_T2)]
EE <- ee[,.(ate,se,lower,upper,A1_T1,A1_T2)]
u <- TT[EE]
u[,.(coverage = mean(lower <= true.ate & true.ate <= upper)),keyby = .(A1_T1,A1_T2)]
u[,mean(ate),by = A1_T2]
u[,table(true.ate)]

######################################################################
### test.R ends here
