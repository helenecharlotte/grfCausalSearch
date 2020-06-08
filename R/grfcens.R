### grfcens.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2020 (16:37) 
## Version: 
## Last-Updated: Jun  8 2020 (13:47) 
##           By: Thomas Alexander Gerds
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
grfcens <- function(formula,
                    data, 
                    CR.as.censoring=FALSE,
                    cause=1,
                    method.weight="KM",
                    formula.weight,
                    args.weight=NULL,
                    treatment,
                    times,...){
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       stripSpecials="intervene",
                                       specials="intervene")
    cens.code <- attr(EHF$event.history,"cens.code")
    dt <- data.table(cbind(unclass(EHF$event.history),EHF$design,EHF$intervene))
    dt[,id:=1:.N]
    if (CR.as.censoring) {
        dt[,is.censored:=as.numeric(event!=cause)]
    }else{
        dt[,is.censored:=as.numeric(status==0)]
    }
    ## foreach::foreach (t0=times) %dopar%{
    t0=times[1]
    if (method.weight=="KM"){
        if (missing(formula.weight)) formula.weight <- formula
        reverse.km.fit <- prodlim(formula.weight, data=data, reverse=TRUE)
        ## resultat of predictSurvIndividual is sorted by (time, -status)
        setorder(dt,time,is.censored)
        #-- define real-valued outcome:
        dt[,Y1:=(event==cause)*(time<=t0)/predictSurvIndividual(reverse.km.fit)]
    } else{
        if (missing(formula.weight))
            formula.weight <- formula(paste0("Surv(time,is.censored)~",paste(c(names(EHF$design),names(EHF$intervene)),collapse="+")))
        else
            formula.weight <- update(formula.weight,"Surv(time,is.censored)~.")
        reverse.forest <- do.call("ranger",c(list(formula=formula.weight, data=dt),args.weight))
        Gmat <- stats::predict(reverse.forest,data=dt)$survival
        jtimes <- ranger::timepoints(reverse.forest)
        Gi.minus <- sapply(1:length(dt$time),function(i){
            pos <- prodlim::sindex(jump.times=jtimes,
                                   eval.times=dt$time[[i]],
                                   comp="smaller",
                                   strict=1L)
            c(1,Gmat[i,])[1+pos]
        })
        dt[,weights:=Gi.minus]
        dt[,Y1:=as.numeric((event==cause)*(time<=t0))]
        dt[Y1!=0,Y1:=1/weights]
    }
    if (any(is.infinite(dt$Y1))) browser()
    ## back to original order
    setorder(dt,id)
    #-- run grf:
    grf.A <- causal_forest(X=EHF$design,
                           Y=dt[["Y1"]],
                           W=as.numeric(as.character(EHF$intervene[[1]])),
                           ...)
    out <- average_treatment_effect(grf.A)
    names(out) <- c("grfate","grfate.se")
    grfate.lower=out[[1]]+qnorm(0.025)*out[[2]]
    grfate.upper=out[[1]]+qnorm(0.975)*out[[2]]
    out <- c(out,grfate.lower=grfate.lower,grfate.upper=grfate.upper)
    out
}
library(riskRegression)
library(prodlim)
library(survival)
library(data.table)
library(grf)
library(ranger)
if (FALSE){
    set.seed(8)
    d <- sampleData(200)
    a <- grfcens(Hist(time,event)~intervene(X4)+X6+X7+X8+X9+X1+X2+X3,data=d,times=10)
    b <- grfcens(Hist(time,event)~intervene(X4)+X6+X7+X8+X9+X1+X2+X3,data=d,times=8,method.weight="ranger",num.trees=50)
    rbind(a,b)
}



######################################################################
### grfcens.R ends here
