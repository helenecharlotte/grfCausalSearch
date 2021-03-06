hunt.fun.fast <- function(dt, ## time, delta (0=censoring), A, X
                          X.vars=paste0("X", 1:10),
                          CR.as.censoring=FALSE,
                          stratify.CR=NULL,
                          stratify.A=FALSE,
                          event=1, method.weight="KM",
                          A.vars=paste0("A", 1:10),
                          t0=0.5,
                          num.trees=1000,
                          browse=FALSE) {

    ATE.list <- lapply(A.vars, function(a) {

        #-- reset weights and define id variable:

        dt[, id:=1:.N]
        dt[, wt:=0]

        if (browse) browser()

        if (length(stratify.CR)>0 & stratify.A) {
            stratify.CR.A <- paste0(stratify.CR, "+", a)
        } else if (stratify.A) {
            stratify.CR.A <- a
        } else {
            stratify.CR.A <- stratify.CR
        }

        #-- compute weights with KM:

        if (method.weight=="KM") {
            if (CR.as.censoring) {
                dt[,is.event:=as.numeric(delta==event)]
            }else{
                dt[,is.event:=delta]
            }
            if (length(stratify.CR.A)==0) {
                form <- formula(paste0("Hist(time, is.event)~ 1"))
            } else{
                form <- formula(paste0("Hist(time, is.event)~ ", stratify.CR.A))
            }
            reverse.km.fit <- prodlim(form, data=dt, reverse=TRUE)
            ## resultat of predictSurvIndividual is sorted by (time, -status)
            setorder(dt,time,-is.event)
            dt[,wt:=(delta==event)*(time<=t0)/predictSurvIndividual(reverse.km.fit)]
            ## back to original order
            setorder(dt,id)
        } else {
            warning("only KM implemented")
        }

        #-- define real-valued outcome:
        dt[, Y1:=wt*as.numeric(time<=t0)]
        #-- run grf:
        
        grf.A <- causal_forest(
            model.matrix(~ -1 + .,
                         dt[, (1:ncol(dt))[colnames(dt) %in% X.vars], with=FALSE]),
            dt[, Y1],
            dt[, (1:ncol(dt))[colnames(dt)==a], with=FALSE][[1]],
            num.trees=num.trees)

        return(average_treatment_effect(grf.A))

    })
    ATE.out <- as.data.frame(do.call("rbind", ATE.list))

    ATE.out$CI.lwr <- ATE.out$estimate - 1.96*ATE.out$std.err
    ATE.out$CI.upr <- ATE.out$estimate + 1.96*ATE.out$std.err

    rownames(ATE.out) <- gsub("A\\.", "", A.vars)

    ATE.out <- ATE.out[order(ATE.out$estimate), ]

    return(ATE.out)
}
