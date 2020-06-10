### grfcens.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2020 (16:37) 
## Version: 
## Last-Updated: Jun 10 2020 (14:53) 
##           By: Thomas Alexander Gerds
##     Update #: 45
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' This function implements a two-step estimator of the average treatment effect 
##' for time-to-event outcome defined as a t-year risk difference
##' according to an intervention. In the presence of competing risks, the function
##' estimates either net or crude probabilities.
##' 
##' One should check the causal three: consistency, positivity and exchangeability.
##' Exchangeability may be a more serious concern when analysing the net probabilities.
##' @title Causal forests for censored data
##' @param formula Formula where left hand side specifies the event history variables
##' and the right hand side the predictor variables where exactly one variable is
##' marked as the exposure of interest via intervene(variable). See examples.
##' @param data Data in which to evaluate the formulae.
##' @param CR.as.censoring Only for settings with competing risks. Logical.
##' If \code{TRUE} analyse net probabilities.
##' @param cause Only for settings with competing risks. The cause of interest.
##' @param method.weight Method for estimating censoring weights.
##' Either \code{'km'} or \code{'ranger'}.
##' @param formula.weight Optional formula specifying the variables used to
##' estimate the censoring weights.
##' @param args.weight Additional arguments passed to method.weight
##' for estimating the censoring weigts.
##' @param times Time point(s) for evaluating the t-year risks.
##' @param ... Additional arguments passed to \code{causal_forest}.
##' @return Average treatment effects
##' @seealso causal_forest
##' @examples
##' library(lava)
##' library(riskRegression)
##' library(prodlim)
##' library(riskRegression)
##' library(prodlim)
##' library(survival)
##' library(data.table)
##' library(grf)
##' library(ranger)
##' m <- lvm(~A+X1+X2+X3+X4+X5)
##' # treatment variable
##' distribution(m,~A) <- binomial.lvm()
##' # latent outcome variables
##' distribution(m,~latenttime1) <- coxWeibull.lvm()
##' distribution(m,~latenttime2) <- coxWeibull.lvm()
##' distribution(m,~censtime) <- coxWeibull.lvm()
##' # observed outcome
##' m <- eventTime(m,time~min(latenttime1=1,latenttime2=2,censtime=0),"event")
##' # dependencies: log-odds ratios and log-hazard ratios
##' # no effect of A on latenttime1
##' regression(m,A~X1+X2+X3+X4+X5) <- c(.1,-.3,.8,-.1,0)
##' regression(m,latenttime1~A+X1+X2+X3+X4+X5) <- c(0,.1,-.3,.8,-.1,0)
##' regression(m,latenttime2~A+X1+X2+X3+X4+X5) <- c(0.2,0,-.3,-.8,-.1,.5)
##' # no effect of any variable on censoring distribution
##' gammaA <- gamma1 <- gamma2 <- gamma3 <- gamma4 <- gamma5  <- 0
##' regression(m,censtime~A+X1+X2+X3+X4+X5) <- c(gammaA,gamma1,gamma2,gamma3,gamma4,gamma5)
##' plot(m)
##' set.seed(99)
##' d <- sim(m,1000)
##' plot(prodlim(Hist(time,event)~A,data=d))
##' fit0 <- prodlim(Hist(time,event)~A,data=d)
##' fit1 <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7)
##' fit1
##' fit2 <- grfcens(formula=Hist(time,event)~intervene(A)+X1+X2+X3+X4+X5,data=d,times=7,num.tree=50,args.weight=list(num.tree=50))
##' fit2
##' @export 
##' @author Helene Charlotte Rytgaard <hely@@sund.ku.dk>, Thomas A. Gerds <tag@@biostat.ku.dk>
grfcens <- function(formula,
                    data, 
                    CR.as.censoring=FALSE,
                    cause=1,
                    method.weight="ranger",
                    formula.weight,
                    args.weight=NULL,
                    treatment,
                    times,...){
    EHF <- prodlim::EventHistory.frame(formula=formula,
                                       data=data,
                                       stripSpecials="intervene",
                                       specials="intervene")
    dt <- data.table(cbind(unclass(EHF$event.history),EHF$design,EHF$intervene))
    if (length(times)>1) warning("Works currently only for one time point at a time.")
    t0=times[[1]]
    cens.code <- attr(EHF$event.history,"cens.code")
    if (attr(EHF$event.history,"cens.type")=="uncensored"){
        dt[,Y1:=as.numeric(time<=t0)]
        if (length(unique(dt[["Y1"]]))<=1) stop(paste0("Outcome status has no variation at time ",t0))
        grf.A <- causal_forest(X=EHF$design,
                               Y=dt[["Y1"]],
                               W=as.numeric(as.character(EHF$intervene[[1]])),
                               ...)
    }else{ ## require right censored
        stopifnot(attr(EHF$event.history,"cens.type")=="rightCensored")
        dt[,id:=1:.N]
        if (CR.as.censoring) {
            dt[,is.censored:=as.numeric(event!=cause)]
        }else{
            dt[,is.censored:=as.numeric(status==0)]
        }
        ## foreach::foreach (t0=times) %dopar%{
        if (tolower(method.weight)=="km"){
            if (missing(formula.weight)){
                formula.weight <- formula(paste0("Hist(time,status)~",paste(c(colnames(EHF$design),colnames(EHF$intervene)),collapse="+")))
                reverse.km.fit <- prodlim(formula.weight, data=dt, reverse=TRUE)
            }else{
                reverse.km.fit <- prodlim(formula.weight, data=data, reverse=TRUE)
            }
            ## resultat of predictSurvIndividual is sorted by (time, -status)
            setorder(dt,time,is.censored)
            #-- define real-valued outcome:
            dt[,Y1:=(event==cause)*(time<=t0)/predictSurvIndividual(reverse.km.fit)]
            if (length(unique(dt[["Y1"]]))<=1) stop(paste0("Outcome status has no variation at time ",t0))
        } else{
            if (missing(formula.weight))
                formula.weight <- formula(paste0("Surv(time,is.censored)~",paste(c(colnames(EHF$design),colnames(EHF$intervene)),collapse="+")))
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
            if (length(unique(dt[["Y1"]]))<=1) stop(paste0("Outcome status has no variation at time ",t0))
        }
        if (any(is.infinite(dt$Y1))) stop(paste0("Weighted outcome status has infinite values at time ",t0))
        ## back to original order
        setorder(dt,id)
        #-- run grf:
        grf.A <- causal_forest(X=EHF$design,
                               Y=dt[["Y1"]],
                               W=as.numeric(as.character(EHF$intervene[[1]])),
                               ...)
    }
    out <- average_treatment_effect(grf.A)
    names(out) <- c("ate","se")
    lower=out[[1]]+qnorm(0.025)*out[[2]]
    upper=out[[1]]+qnorm(0.975)*out[[2]]
    out <- c(out,lower=lower,upper=upper)
    out
}

######################################################################
### grfcens.R ends here
