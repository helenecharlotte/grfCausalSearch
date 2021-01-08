### causalhunter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2020 (16:37) 
## Version: 
## Last-Updated: Jun 27 2020 (09:15) 
##           By: Thomas Alexander Gerds
##     Update #: 117
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
##' mydata <- sampleData(200,outcome="competing.risk")
##' weighter(Hist(time,event)~X1+X5+X7,data=mydata,times=5)
##' causalhunter(Hist(time,event)~intervene(X1)+X5+X7,data=mydata,times=5)
##' causalhunter(Hist(time,event)~intervene(X1)+intervene(X2)+X5+X7,data=mydata,times=5)
##' causalhunter(Hist(time,event)~intervene(X1)+intervene(X2)+X5+X7,data=mydata,times=c(5,8))
##' 
##' @export 
##' @author Helene Charlotte Rytgaard <hely@@sund.ku.dk>, Thomas A. Gerds <tag@@biostat.ku.dk>
causalhunter <- function(formula,
                         data, 
                         CR.as.censoring=FALSE,
                         cause=1,
                         method.weight="ranger",
                         fit.separate=FALSE,
                         formula.weight,
                         args.weight=NULL,
                         times,truncate=TRUE,...){
    EHF <- prodlim::EventHistory.frame(formula=formula,
                                       data=data,
                                       stripSpecials="intervene",
                                       specials="intervene",
                                       specialsFactor=FALSE,
                                       unspecialsDesign=TRUE) # grf does not work with non-numeric features
    if (missing(formula.weight) || is.null(formula.weight)){
        variables <- unique(c(c(names(attr(EHF$design, "levels")),
                                colnames(EHF$design))[c(names(attr(EHF$design, "levels")),
                                                        colnames(EHF$design))%in%colnames(data)],
                              colnames(EHF$intervene)))
        formula.weight <- formula(paste0("Hist(time,event)~",paste(variables,collapse="+")))
        dt <- data.table(cbind(unclass(EHF$event.history),
                               #EHF$design
                               data[, variables, with=FALSE],EHF$intervene))
        ## fix censored event status 
        dt[status==0,event:=0]
        Y <- do.call("weighter",c(list(formula=formula.weight,data=dt,times=times,
                                       method=method.weight,CR.as.censoring=CR.as.censoring,
                                       fit.separate=fit.separate,
                                       cause=cause,truncate=truncate),args.weight))
    }else{
        Y <- do.call("weighter",c(list(formula=formula.weight,data=data,times=times,
                                       method=method.weight,CR.as.censoring=CR.as.censoring,
                                       fit.separate=fit.separate,
                                       cause=cause,truncate=truncate),args.weight))
    }
    if (any(is.infinite(Y))) stop(paste0("Weighted outcome status has infinite values at time ",t0))
    if (length(unique(Y))<=1) stop(paste0("Outcome status has no variation at time ",t0))
    #-- run grf:
    out <- NULL
    for (s in 1:length(times)){
        for (i in 1:NCOL(EHF$intervene)){
            grf.A <- causal_forest(X=EHF$design,
                                   Y=Y[,s,drop=TRUE],
                                   W=as.numeric(as.character(EHF$intervene[[i]])),
                                   ...)
            est <- average_treatment_effect(grf.A)
            names(est) <- c("ate","se")
            lower=est[[1]]+qnorm(0.025)*est[[2]]
            upper=est[[1]]+qnorm(0.975)*est[[2]]
            est <- data.frame(time=times[s],intervene=colnames(EHF$intervene)[[i]],ate=est[[1]],se=est[[2]],lower=lower,upper=upper)
            out <- rbind(out,est)
        }
    }
    rownames(out) <- NULL
    out
}

######################################################################
### causalhunter.R ends here
