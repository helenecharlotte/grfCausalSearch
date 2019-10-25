plot.output <- function(vimp.obj, title=bquote(hat(theta)[1] ~ "(unadjusted weights)"), which=1,
                        sim=FALSE, text.pos=0.01, text.pos.cov=0.01, psi=1) {

    if (sim) {

        psi0 <- vimp.obj[[1]]

        coverage <- unlist(lapply(1:n.A, function(a) {
            return(mean(unlist(lapply(vimp.obj[-1], function(est) {
                est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
                est.a <- est[a, ]
                return(as.numeric(est.a$CI.lwr<=psi0[a] & psi0[a]<=est.a$CI.upr))
            }))))
        }))

        mean <- unlist(lapply(1:n.A, function(a) {
            return(mean(unlist(lapply(vimp.obj[-1], function(est) {
                est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
                est.a <- est[a, ]
                return(est.a$estimate)
            }))))
        }))

        CI.lwr <- unlist(lapply(1:n.A, function(a) {
            return(mean(unlist(lapply(vimp.obj[-1], function(est) {
                est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
                est.a <- est[a, ]
                return(est.a$CI.lwr)
            }))))
        }))

        CI.upr <- unlist(lapply(1:n.A, function(a) {
            return(mean(unlist(lapply(vimp.obj[-1], function(est) {
                est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
                est.a <- est[a, ]
                return(est.a$CI.upr)
            }))))
        }))

        sd <- unlist(lapply(1:n.A, function(a) {
            return(sd(unlist(lapply(vimp.obj[-1], function(est) {
                est <- est[order(as.numeric(gsub("A", "", rownames(est)))),]
                est.a <- est[a, ]
                return(est.a$estimate)
            }))))
        }))
        
        
        dt.vimp <- data.table(treatment=rownames(vimp.obj[[2]][order(as.numeric(
                                                             gsub("A", "", rownames(vimp.obj[[2]])))),]),
                              estimate=mean,
                              CI.lwr=mean-1.96*sd,
                              CI.upr=mean+1.96*sd,
                              psi0=psi0,
                              coverage=coverage)

        
    } else {

        vimp.obj.sorted <- #sort.fun(
            vimp.obj#)
        
        dt.vimp <- data.table(treatment=rownames(vimp.obj.sorted),
                              vimp.obj.sorted)

    }


    dt.vimp[, y.var:=1:.N]
    dt.vimp[, col:=factor(1*(CI.lwr<0 & CI.upr>0)+2*(CI.upr<0))]

    if (sim) {
        
        dt.vimp.psi <- dt.vimp[psi,]
        dt.zero <- data.table(x=c(0, dt.vimp.psi[, psi0]),
                              y=rep(-0.5, length=1+nrow(dt.vimp.psi)),
                              y.var=max(dt.vimp[, y.var]),
                              label=c("0", paste0("bar(theta)[", which, "*','*A[", psi,"]]", sep="")))

    } else {

        dt.zero <- data.table(x=0, y=-0.5, label="0")
        
    }



    colors.manual <- 1:3
    if (length(unique(dt.vimp[, col]))==2) colors.manual <- c(2,3)

    p.out <- ggplot(dt.vimp) + theme_void() +
        geom_point(aes(x=estimate, y=y.var, col=col), size=2) +
        geom_segment(aes(x=min(CI.lwr)-0.003, xend=CI.lwr, y=y.var, yend=y.var), col="gray93", size=0.1) +
        geom_segment(aes(x=CI.lwr, xend=CI.upr,
                         y=y.var, yend=y.var, col=col), size=0.7) +
        geom_segment(aes(x=CI.lwr, xend=CI.lwr,
                         y=y.var-0.15, yend=y.var+0.15, col=col), size=0.7) +
        geom_segment(aes(x=CI.upr, xend=CI.upr,
                         y=y.var-0.15, yend=y.var+0.15, col=col), size=0.7) +
        geom_text(aes(x=min(CI.lwr)-text.pos, y=y.var, label=treatment)) +
        geom_segment(data=dt.zero, aes(x=x, xend=x, y=0, yend=y.var+1), linetype=2) +
        geom_text(data=dt.zero, aes(x=x, y=y, label=label), parse=TRUE) +
        geom_segment(aes(x=min(CI.lwr), xend=max(CI.upr), y=0, yend=0)) +
        geom_segment(aes(x=min(CI.lwr)-0.006, xend=min(CI.lwr)-0.006, y=1-0.7, yend=max(y.var)+0.5)) +
        geom_segment(aes(x=min(CI.lwr)-0.006, xend=min(CI.lwr)-0.007, y=1-0.7, yend=1-0.7)) +
        geom_segment(aes(x=min(CI.lwr)-0.006, xend=min(CI.lwr)-0.007, y=max(y.var)+0.5, yend=max(y.var)+0.5)) +
        geom_segment(aes(x=min(CI.lwr), xend=min(CI.lwr), y=0-0.24, yend=0)) +
        geom_segment(aes(x=max(CI.upr), xend=max(CI.upr), y=0-0.24, yend=0)) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.x=element_blank(),
              plot.background=element_blank(),
              legend.position="none",
              plot.title=element_text(size=14, hjust=0.5)) +
        scale_color_manual(values=hue_pal()(3)[colors.manual]) +
        ggtitle(title)

    if (sim) {

        p.out <- p.out +
            geom_text(aes(x=max(CI.upr)+text.pos.cov, y=y.var, label=paste0(formatC(round(coverage*100,1), format='f', digits=1), " %")))
            
    }

    return(p.out)
}


sort.fun <- function(ATE.out) {
    return( ATE.out[order(ATE.out$estimate), ])
}
