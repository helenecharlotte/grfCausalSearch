### weighter.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun 27 2020 (06:33) 
## Version: 
## Last-Updated: Jan 11 2021 (19:09) 
##           By: Thomas Alexander Gerds
##     Update #: 63
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
                     times,
                     truncate=FALSE,
                     ...){
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
            reverse.km.fit <- prodlim(ff, data=dt, reverse=TRUE)
            ## result of predictSurvIndividual is sorted by (time, -status)
            dt[,id:=1:.N]
            setorder(dt,time,is.censored)
            dt[,Weight:=(event==cause)/predictSurvIndividual(reverse.km.fit)]
            ## back to original order
            setorder(dt,id)
            Weight <- dt[["Weight"]]
        } else {
            if (fit.separate & CR.as.censoring) {
                dt[,Weight:=1]
                for (jj in c(ifelse(CR.as.censoring, c(0,2), 0))) {
                    ff <- update(formula,paste0("Surv(time,event==", jj, ")~."))
                    reverse.forest <- do.call("ranger",c(list(formula=ff, data=dt),...))
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
                args <- c(list(formula=ff, data=dt),list(...))
                N <- NROW(dt)
                tune.node.size <- length(args$min.node.size)
                if (tune.node.size>1){
                    tune.node.size <- args$min.node.size
                    stopifnot(all(tune.node.size<N))
                    reverse.forest.candidates <- lapply(c(tune.node.size),function(nodesize){
                        rargs <- args
                        rargs$min.node.size <- nodesize
                        rargs$tuning.time <- NULL
                        do.call("ranger",rargs)})
                    names(reverse.forest.candidates) <- paste0("MNS:",tune.node.size)
                    suppressMessages(x <- Score(reverse.forest.candidates,times=args$tuning.time,formula=ff,metrics="brier",data=dt,split.method="cv10",se.fit=FALSE,contrast=NULL,progress.bar=NULL))
                    winner <- as.character(setkey(x$Brier$score,Brier)[1][["model"]])
                }else{
                    # no tuning of nodesize
                    reverse.forest.candidates <- lapply(1,function(nodesize){
                        rargs <- args
                        rargs$tuning.time <- NULL
                        do.call("ranger",rargs)
                    })
                    names(reverse.forest.candidates) <- paste0("MNS:",args$min.node.size)
                    winner <- names(reverse.forest.candidates)
                }
                # print(winner)
                if (winner=="Null model"){ # marginal Kaplan-Meier
                    if (CR.as.censoring==TRUE) dt[is.censored==1,status:=0]
                    reverse.km.fit <- prodlim(Hist(time,status)~1, data=dt, reverse=TRUE)
                    ## result of predictSurvIndividual is sorted by (time, -status)
                    dt[,id:=1:.N]
                    setorder(dt,time,is.censored)
                    dt[,Weight:=(event==cause)/predictSurvIndividual(reverse.km.fit)]
                    ## back to original order
                    setorder(dt,id)
                    Weight <- dt[["Weight"]]
                }else{
                    reverse.forest <- reverse.forest.candidates[[winner]]
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
    }
    Y <- do.call("cbind",lapply(times,function(t0){as.numeric(dt[["time"]]<=t0)*Weight}))
    if (method=="ranger") attr(Y,"winner") <- winner else attr(Y,method)
    Y
}


######################################################################
### weighter.R ends here
