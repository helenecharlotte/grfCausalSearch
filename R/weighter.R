### weighter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 27 2020 (06:33) 
## Version: 
## Last-Updated: Jun 27 2020 (10:09) 
##           By: Thomas Alexander Gerds
##     Update #: 19
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
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=TRUE)
##' weighter(formula=Hist(time,event)~X1+X5+X7,data=mydata,times=c(5,8),truncate=TRUE,method="km")
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
                     method="ranger",
                     times,truncate=FALSE,...){
    EHF <- prodlim::EventHistory.frame(formula=formula,
                                       data=data,
                                       specials=NULL,
                                       unspecialsDesign=FALSE)
    response.type <- attr(EHF$event.history,"model") # either "survival" or "competing.risk"
    cens.type <- attr(EHF$event.history,"cens.type") # either "uncensored" or "rightCensored"
    if (CR.as.censoring == TRUE && response.type=="survival")
        warning("Only one cause of event in data. Removing a competing risk which does not occur.")
    dt <- data.table(cbind(unclass(EHF$event.history),EHF$design))
    ## set censored to 0
    dt[,event:=event*status]
    if (truncate[[1]]==TRUE){
        dt[time>max(times),event:=0]
        dt[time>max(times),time:=max(times)]
    }
    # uncensored
    if (CR.as.censoring==FALSE && cens.type=="uncensored"){
        Weight <- rep(1,NROW(data))
    }else{
        ## require either right censored data, or request for the hypothetical world
        ## where competing risk has been removed
        stopifnot(cens.type=="rightCensored" || CR.as.censoring == TRUE)
        ## in case without competing risks add event variable
        if (response.type=="survival") dt[,event:=status]
        if (CR.as.censoring) {
            dt[,is.censored:=as.numeric(event!=cause)]
        }else{
            dt[,is.censored:=as.numeric(status==0)]
        }
        if (tolower(method)=="km"){
            if (CR.as.censoring==TRUE) dt[is.censored==1,status:=0]
            ff <- update(formula,"Hist(time,status)~.")
            if (FALSE) {
                coef(glmnet(y=dt[, Surv(time,event!=1)], x=as.matrix(EHF$design),
                            family="cox", maxit=1000, lambda=seq(0.0001, 0.1, length=10)))
                coef(glmnet(y=dt[, Surv(time,event==2)], x=as.matrix(EHF$design),
                            family="cox", maxit=1000, lambda=seq(0.0001, 0.1, length=10)))
            }
            reverse.km.fit <- prodlim(ff, data=dt, reverse=TRUE)
            ## result of predictSurvIndividual is sorted by (time, -status)
            dt[,id:=1:.N]
            setorder(dt,time,is.censored)
            dt[,Weight:=(event==cause)/predictSurvIndividual(reverse.km.fit)]
            ## back to original order
            setorder(dt,id)
            Weight <- dt[["Weight"]]
        } else{
            #ff <- update(formula,"Surv(time,is.censored)~.")
            if (fit.separate & CR.as.censoring) {
                dt[,Weight:=1]
                for (jj in c(ifelse(CR.as.censoring, c(0,2), 0))) {
                    ff <- update(formula,paste0("Surv(time,event==", jj, ")~."))
                    reverse.forest <- do.call("ranger",c(list(formula=ff, data=dt),...))
                    if (FALSE) {
                        ff <- update(formula,paste0("Surv(time,event==", 2, ")~."))
                        test <- do.call("ranger",c(list(formula=ff, data=dt, importance="permutation"),...))
                        sort(importance(test))
                    }
                    Gmat <- stats::predict(reverse.forest,data=dt)$survival
                    jtimes <- ranger::timepoints(reverse.forest)
                    Gi.minus <- sapply(1:length(dt$time),function(i){
                        pos <- prodlim::sindex(jump.times=jtimes,eval.times=dt$time[[i]],
                                               comp="smaller",strict=1L)
                        c(1,Gmat[i,])[1+pos]
                    })
                    dt[,Weight:=1/Gi.minus]
                    Weight <- dt[["Weight"]]
                }
                dt[event!=cause,Weight:=0]
                
            } else {
                ff <- update(formula,"Surv(time,is.censored)~.")
                reverse.forest <- do.call("ranger",c(list(formula=ff, data=dt),...))
                Gmat <- stats::predict(reverse.forest,data=dt)$survival
                jtimes <- ranger::timepoints(reverse.forest)
                Gi.minus <- sapply(1:length(dt$time),function(i){
                    pos <- prodlim::sindex(jump.times=jtimes,eval.times=dt$time[[i]],comp="smaller",strict=1L)
                    c(1,Gmat[i,])[1+pos]
                })
                dt[,Weight:=1/Gi.minus]
                dt[event!=cause,Weight:=0]
                Weight <- dt[["Weight"]]
            }
        }
    }
    Y <- do.call("cbind",lapply(times,function(t0){as.numeric(dt[["time"]]<=t0)*Weight}))
    Y
}


######################################################################
### weighter.R ends here
