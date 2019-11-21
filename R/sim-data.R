sim.data <- function(n=1000,
                     compute.counterfactuals=TRUE,
                     compute.psi0=c(0,1,2),
                     which.A=1,
                     seed=runif(1, min=1, max=30000),
                     CR=c(1,2),
                     C.shape=0.4,
                     t0=0.5, n.A=10, cens.A=0,
                     form.A = function(X, alpha=0.2, beta=0.5) 0.5+alpha+as.numeric(X[,1])*beta,
                     form.T1 = function(X, A) -1.1 + as.numeric(X[, 1])*0.2 -
                         as.numeric(X[,3])*0.1 -
                         A[, 1]*1.5,
                     form.T2 = function(X, A) 0.1 - as.numeric(X[, 2])*0.4 -
                                              as.numeric(X[, 1])*0.33 +
                                              1.5*A[, 2]) {

    if (length(seed)>0) set.seed(seed)

    rexpit <-  function(x) rbinom(n=length(x), size=1, prob=plogis(x))

    X <- data.frame(X1 = runif(n),
                    X2 = factor(sample(1:3, n, TRUE)),
                    X3 = factor(sample(1:4, n, TRUE)),
                    X4 = runif(n),
                    X5 = runif(n),
                    X6 = runif(n))

    A <- data.frame(sapply(1:n.A, function(a) {
        return(rexpit(form.A(X[, rep((a%%3)*1+(a%%2)*1+1-1*(a==5), 4)], alpha=(2*(a%%2)-1)/(n.A*10),
                             beta=((2*(a%%2)-1)+((a%%3)>0))/10)))
    }))

    names(A) <- paste0("A", 1:n.A)

    if (cens.A>0) {
        C <- rweibull(n, shape=C.shape+A[, cens.A], scale=2.4)
    } else {
        C <- rweibull(n, shape=C.shape+0.5, scale=2.4)
    }

    if (compute.counterfactuals) {

        A.1 <- A.0 <- A
        A.1[, which.A] <- 1
        A.0[, which.A] <- 0

        T1.1 <- rweibull(n, shape=0.8-0.5*form.T1(X, A.1), scale=1)
        T1.0 <- rweibull(n, shape=0.8-0.5*form.T1(X, A.0), scale=1)
        T1   <- T1.1*A[, which.A] + T1.0*(1-A[, which.A])

        if (CR[1]>1) {

            T2.1 <- rweibull(n, shape=0.8-0.5*form.T2(X, A.1), scale=1)
            T2.0 <- rweibull(n, shape=0.8-0.5*form.T2(X, A.0), scale=1)
            T2   <- T2.1*A[, which.A] + T2.0*(1-A[, which.A])

            T.1  <- sapply(1:n, function(i) min(T1.1[i], T2.1[i]))
            T.0  <- sapply(1:n, function(i) min(T1.0[i], T2.0[i]))
            time <- sapply(1:n, function(i) min(T1[i], T2[i], C[i]))

            delta.1 <- 1*(T.1==T1.1) + 2*(T.1==T2.1)
            delta.0 <- 1*(T.0==T1.0) + 2*(T.0==T2.0)
            delta   <- 1*(time==T1) + 2*(time==T2)

            df <- cbind(X, A, time=time, delta=delta,
                        T.1=T.1, delta.1, T.0=T.0, delta.0,
                        T1.1=T1.1, T1.0=T1.0, T2.1=T2.1, T2.0=T2.0)

        } else {

            time  <- sapply(1:n, function(i) min(T1[i], C[i]))
            delta <- 1*(time==T1)

            df <- cbind(X, A, time=time, delta=delta,
                        T1=T1, T1.1=T1.1, T1.0=T1.0)

        }
    } else {

        T1 <- rweibull(n, shape=0.8-0.5*form.T1(X, A), scale=1)

        if (CR[1]>1) {

            T2    <- rweibull(n, shape=0.8-0.5*form.T2(X, A.1), scale=1)
            time  <- sapply(1:n, function(i) min(T1[i], T2[i], C[i]))
            delta <- 1*(time=T1) + 2*(time=T2)
            df    <- cbind(X, A, time=time, delta=delta, T1=T1, T2=T2)

        } else {

            time  <- sapply(1:n, function(i) min(T1[i], C[i]))
            delta <- 1*(time=T1)
            df <- cbind(X, A, time=time, delta=delta, T1=T1)

        }
    }

    if (compute.psi0[1]>0 & CR[1]>1) {
        if (compute.psi0[1]==1) {
            return(mean((df$T1.1<=t0)-(df$T1.0<=t0)))
        } else {
            return(mean((df$delta.1==1 & df$T.1<=t0)-(df$delta.0==1 & df$T.0<=t0)))
        }
    } else if (compute.psi0[1]>0 & CR[1]==1) {
        return(mean(df$T1.1<=t0) - mean(df$T1.0<=t0))
    } else {
        return(setDT(df))
    }

}
