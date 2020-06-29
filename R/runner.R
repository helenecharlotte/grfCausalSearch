### runner.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 27 2020 (09:36) 
## Version: 
## Last-Updated: Jun 29 2020 (08:26) 
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
runner <- function(seed,
                   cens=0.2,
                   effect.A2=1.5,
                   M=200,
                   n=200,
                   NT=50,
                   args.weight=list(num.tree=50,replace=FALSE,probability=TRUE),
                   method.weight="ranger",
                   formula.weight=NULL,
                   truncate=TRUE,
                   intervene="all",
                   effect="net",
                   verbose=TRUE,
                   cores=25){
    if (cores>1) registerDoParallel(cores=cores) else registerDoSEQ()
    if (verbose) pb <- txtProgressBar(max = M, style = 3,width=20)
    v <- seed
    ## all net-effects are zero except for A1
    truth <- data.table(intervene=paste0("A",1:10),truth=0)
    # mean(sapply(1:10,function(u){sim.data(n=1000000,seed=u,compute.psi=1,CR=2,C.shape=1,which.A=1,form.T1 = function(X, A){ -1.1 + as.numeric(X[, 1])*0.2 - as.numeric(X[,3])*0.1 - A[, 1]*1.5},form.T2 = function(X, A) {0.1 - as.numeric(X[, 2])*0.4 - as.numeric(X[, 1])*0.33 + 1.5*A[, 2]})}))
    truth[1,truth:=-0.1122533]
    #truth[1,truth:=sim.data(n=1000000,seed=seed,compute.psi=1,CR=2,C.shape=cens,which.A=1,form.T1 = function(X, A){ -1.1 + as.numeric(X[, 1])*0.2 - as.numeric(X[,3])*0.1 - A[, 1]*1.5},form.T2 = function(X, A) {0.1 - as.numeric(X[, 2])*0.4 - as.numeric(X[, 1])*0.33 + effect.A2*A[, 2]})]
    ## print(truth)
    result <- foreach(s=1:M,.combine="rbind",.errorhandling="pass")%dopar%{
        if (verbose) setTxtProgressBar(pb, s)
        d <- sim.data(n=n,seed=seed+s,compute.psi=0,CR=2,C.scale=cens,which.A=2,form.T1 = function(X, A){ -1.1 + as.numeric(X[, 1])*0.2 - as.numeric(X[,3])*0.1 - A[, 1]*1.5},form.T2 = function(X, A) {0.1 - as.numeric(X[, 2])*0.4 - as.numeric(X[, 1])*0.33 + effect.A2*A[, 2]})
        set.seed(seed+s)
        if (intervene=="all"){
            ff <- Hist(time,delta)~intervene(A1)+intervene(A2)+intervene(A3)+intervene(A4)+intervene(A5)+intervene(A6)+intervene(A7)+intervene(A8)+intervene(A9)+intervene(A10)+X1+X2+X3+X4+X5+X6
        } else{
            if (intervene=="A2"){
                ff <- Hist(time,delta)~A1+intervene(A2)+X1+X2+X3+X4+X5+X6
            } else{ 
                ff <- Hist(time,delta)~intervene(A1)+intervene(A2)+X1+X2+X3+X4+X5+X6
            }
        }
        E <- causalhunter(formula=ff,
                          method.weight=method.weight,
                          formula.weight=formula.weight,
                          CR.as.censoring=(effect=="net"),
                          num.tree=NT,
                          args.weight=args.weight,
                          truncate=truncate,
                          data=d,
                          times=0.5)
        if (intervene=="all") E <- cbind(E,rank=rank(-abs(E[["ate"]])))
        cbind(cens=mean(d$time<=0.5&d$delta==0),E)
    }
    if (cores>1) registerDoSEQ()
    cat("\n")
    result <- data.frame(result)
    setDT(result)
    setkey(result,intervene)
    setkey(truth,intervene)
    result <- truth[result]
    result[,coverage:=(lower<truth)&(upper>truth)]
    if (intervene=="all")
        result[,.(time=time[1],
                  cens=mean(cens,na.rm=TRUE),
                  truth=truth[1],
                  num.trees=NT,
                  mean=mean(ate,na.rm=TRUE),
                  se=sd(ate,na.rm=TRUE),
                  mean.se=mean(se,na.rm=TRUE),
                  bias=mean(truth-ate,na.rm=TRUE),
                  abs.bias=mean(abs(truth-ate),na.rm=TRUE),
                  coverage=mean(coverage,na.rm=TRUE),
                  mean.rank=mean(rank,na.rm=TRUE)),by=intervene]
    else
        result[,.(time=time[1],
                  cens=mean(cens,na.rm=TRUE),
                  truth=truth[1],
                  num.trees=NT,
                  mean=mean(ate,na.rm=TRUE),
                  se=sd(ate,na.rm=TRUE),
                  mean.se=mean(se,na.rm=TRUE),
                  bias=mean(truth-ate,na.rm=TRUE),
                  abs.bias=mean(abs(truth-ate),na.rm=TRUE),
                  coverage=mean(coverage,na.rm=TRUE)),by=intervene]
}

## x <- runner(seed=7,cens=.2,effect.A2=1.5,M=25,cores=25,intervene="A2")

## x <- runner(seed=7,cens=.2,effect.A2=1.5,M=2,cores=1,intervene="all")

## x <- runner(seed=7,cens=.2,effect.A2=1.5,M=25,cores=25,intervene="all")


## x <- runner(seed=7,cens=.2,effect.A2=1.5,M=1000,cores=25)
## x
    ## intervene time     cens     truth num.trees          mean         se    mean.se          bias   abs.bias coverage
 ## 1:        A1  0.5 0.199325 -0.112243        50 -0.0976049928 0.08987106 0.09230855 -0.0146380072 0.07268399    0.928
 ## 2:        A2  0.5 0.199325  0.000000        50 -0.0475179808 0.07193446 0.07733154  0.0475179808 0.06837347    0.927
 ## 3:        A3  0.5 0.199325  0.000000        50  0.0011735517 0.07313254 0.07871713 -0.0011735517 0.05792535    0.958
 ## 4:        A4  0.5 0.199325  0.000000        50 -0.0033012148 0.07336529 0.07672876  0.0033012148 0.05912690    0.962
 ## 5:        A5  0.5 0.199325  0.000000        50 -0.0011834901 0.08153129 0.08446799  0.0011834901 0.06554827    0.942
 ## 6:        A6  0.5 0.199325  0.000000        50  0.0004121832 0.07202688 0.07526722 -0.0004121832 0.05632123    0.952
 ## 7:        A7  0.5 0.199325  0.000000        50 -0.0032827355 0.08365978 0.08471408  0.0032827355 0.06678488    0.944
 ## 8:        A8  0.5 0.199325  0.000000        50  0.0004297848 0.07393418 0.07609039 -0.0004297848 0.05858358    0.952
 ## 9:        A9  0.5 0.199325  0.000000        50  0.0022473395 0.07978461 0.07869843 -0.0022473395 0.06345135    0.937
## 10:       A10  0.5 0.199325  0.000000        50 -0.0046262086 0.07389959 0.07681628  0.0046262086 0.05875851    0.953

## y <- runner(seed=7,cens=.2,effect.A2=1.5,M=1000,cores=25,formula.weight=Hist(time,delta)~A2+X1+X2)
## y
    ## intervene time     cens     truth num.trees          mean         se    mean.se          bias   abs.bias coverage
 ## 1:        A1  0.5 0.199325 -0.112243        50 -0.1019012705 0.10121179 0.10110555 -0.0103417295 0.08041739    0.922
 ## 2:        A2  0.5 0.199325  0.000000        50 -0.0243176517 0.07497224 0.08254679  0.0243176517 0.06257119    0.959
 ## 3:        A3  0.5 0.199325  0.000000        50  0.0022570534 0.08223133 0.08664070 -0.0022570534 0.06532132    0.954
 ## 4:        A4  0.5 0.199325  0.000000        50 -0.0034371707 0.08360648 0.08478430  0.0034371707 0.06686156    0.955
 ## 5:        A5  0.5 0.199325  0.000000        50 -0.0008476283 0.09211579 0.09290033  0.0008476283 0.07408740    0.939
 ## 6:        A6  0.5 0.199325  0.000000        50  0.0005090935 0.08102814 0.08312656 -0.0005090935 0.06366563    0.952
 ## 7:        A7  0.5 0.199325  0.000000        50 -0.0023011407 0.09334548 0.09306043  0.0023011407 0.07461658    0.939
 ## 8:        A8  0.5 0.199325  0.000000        50 -0.0010154538 0.08401025 0.08417076  0.0010154538 0.06610295    0.955
 ## 9:        A9  0.5 0.199325  0.000000        50  0.0021362394 0.08892307 0.08679779 -0.0021362394 0.07088711    0.932
## 10:       A10  0.5 0.199325  0.000000        50 -0.0030171077 0.08216301 0.08444099  0.0030171077 0.06521619    0.950


######################################################################
### runner.R ends here
