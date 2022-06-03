### weighter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 27 2020 (06:33) 
## Version: 
## Last-Updated: Apr 29 2022 (12:38) 
##           By: Thomas Alexander Gerds
##     Update #: 126
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param formula
##' @param data
##' @param CR.as.censoring
##' @param cause
##' @param method
##' @param times
##' @param ...
##' @return 
##' @seealso 
##' @examples
##' library(riskRegression)
##' library(survival)
##' library(prodlim)
##' library(data.table)
##' library(grf)
##' library(ranger)
##' mydata <- sampleData(200,outcome="competing.risk")
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=5)
##' set.seed(8)
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=FALSE)
##' set.seed(8)
##' u=weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=TRUE)
##' v=weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=TRUE,method="km")
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=FALSE,method="km")
##' a <- weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=5)
##' b <- weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=5,CR.as.censoring=TRUE)
##' cbind(a,b)
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
weighter <- function(formula,
                     data,
                     CR.as.censoring=FALSE,
                     fit.separate=FALSE,
                     cause=1,
                     method="km",
                     times,
                     truncate=FALSE,
                     ...){
    # parse formula in data
    EHF <- prodlim::EventHistory.frame(formula=formula,data=data,specials=NULL,unspecialsDesign=FALSE)
    X <- EHF$design
    response.type <- attr(EHF$event.history,"model") # either "survival" or "competing.risk"
    if(response.type!="competing.risk") {
        CR.as.censoring <- FALSE
        fit.separate <- FALSE
    }
    cens.type <- attr(EHF$event.history,"cens.type") # either "uncensored" or "rightCensored"
    if (CR.as.censoring == TRUE)
        warning("Only one cause of event in data. Removing a competing risk which does not occur!?")
    dt <- data.table(cbind(unclass(EHF$event.history),X))
    ## set censored to have event value 0
    if (response.type=="survival") dt[,event:=status] else dt[,event:=event*status]
    ## truncate all times that are larger than the maximum evaluation time
    if (truncate[[1]]==TRUE){
        dt[time>max(times),event:=0]
        dt[time>max(times),time:=max(times)+max(times)/10000]
    }
    # case: uncensored and not intervening on competing risk
    if (CR.as.censoring==FALSE && cens.type=="uncensored"){
        Weight <- as.numeric(dt[["event"]] == cause)
    }else{
        ## cases: either right censored data or net effects in the hypothetical world
        ## where the competing risk has been removed
        stopifnot(cens.type=="rightCensored" || CR.as.censoring == TRUE)
        if (CR.as.censoring) {
            dt[,is.censored:=as.numeric(event!=cause)]
        }else{
            dt[,is.censored:=as.numeric(status==0)]
        }
        dt[,Weight:=1]
        if (fit.separate) { JJ <- c(0,2)} else {JJ <- 1}
        if (tolower(method)=="km"){
            for (jj in JJ){
                # reverse Kaplan-Meier for censoring distribution 
                if (jj == 1)
                    dt[,Status := as.numeric(is.censored != jj)]
                else
                    dt[,Status := as.numeric(is.censored == jj)]
                ff <- update(formula,"Hist(time,Status)~.")
                reverse.km.fit <- prodlim(ff, data=dt, reverse=TRUE)
                ## result of predictSurvIndividual is sorted by (time, -status)
                dt[,id:=1:.N]
                setorder(dt,time,is.censored)
                dt[,Weight:=Weight*(event==cause)/predictSurvIndividual(reverse.km.fit)]
                ## back to original order
                setorder(dt,id)
            }
            Weight <- dt[["Weight"]]
        } else {
            # method: ranger
            for (jj in JJ) {
                if (jj == 1)
                    ff <- update(formula,paste0("Surv(time,event!=", jj, ")~."))
                else
                    ff <- update(formula,paste0("Surv(time,event==", jj, ")~."))
                reverse.forest <- do.call("ranger",c(list(formula=ff,
                                                          data=dt,
                                                          replace = FALSE,
                                                          splitrule = "maxstat",
                                                          alpha = 0.05),...))
                Gmat <- stats::predict(reverse.forest,data=dt)$survival
                jtimes <- ranger::timepoints(reverse.forest)
                Gi.minus <- sapply(1:length(dt$time),function(i){
                    pos <- prodlim::sindex(jump.times=jtimes,eval.times=dt$time[[i]],
                                           comp="smaller",strict=1L)
                    c(1,Gmat[i,])[1+pos]
                })
                # multiply the weights G and G_2 
                dt[,Weight:=Weight*(event == cause)/Gi.minus]
            }
            Weight <- dt[["Weight"]]
        }
    }
    Y <- do.call("cbind",lapply(times,function(t0){
        # weigth is already zero for event!=cause
        y = as.numeric(dt[["time"]]<=t0)*Weight
        if (any(is.infinite(y))) stop(paste0("Weighted outcome status has infinite values at time ",t0))
        if (length(unique(y))<=1) {
            print(table(y))
            print(head(dt[["time"]]))
            print(head(Weight))
            stop(paste0("Outcome status has no variation at time ",t0))
        }
        y
    }))
    Y
}


######################################################################
### weighter.R ends here
