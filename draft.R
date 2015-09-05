library(ggplot2)

# sqrtFitMSE <- function(q, ell, N) {
#         p <- length(ell)
#         x <- ell[(q+1):(q+N)]
#         u <- ((p-q):(p-q-N+1))/(p-q)
#         y <- (1-u)^(2/3)
#         fit <- lm(x~y)
#         mu <- -fit$coefficients[1] / fit$coefficients[2]
#         return(list(mse = mean(fit$residuals^2), mu = mu))
# }
# 
# constrainedSqrtFitMSE <- function(q, ell, N) {
#         p <- length(ell)
#         x <- ell[(q+1):(q+N)]
#         u <- ((p-q):(p-q-N+1))/(p-q)
#         y <- (1-u)^(2/3)
#         l <- ell[q+2]
#         if (q==0) {
#                 L <- Inf 
#         } else {
#                 L <- ell[q]
#         }
#         constrainedFit <- nls(x ~ m + alpha*y, algorithm="port",
#             start=c(m=ell[q+1], alpha=1),
#             lower=c(m=l,alpha=-Inf),
#             upper=c(m=L,alpha=Inf)
#             )
#         return(mean(resid(constrainedFit)^2))
# }

advancedSqrtFitMSE <- function(q, ell, N, d=10) {
        p <- length(ell)
        i <- p - which.max(ell[p:1]>0) + 1 #index of the last non-zero value at ell
        range <- (ell[q+1] - ell[i]) / d
        j <- which.min(ell > (ell[q+1]-range))
        if (j < q+10) {j <- q+10}
        x <- ell[(q+1):j]
        #        u <- ((p-q):(p-j+1))/(p-q)
        u <- ((p-q):(p-j+1))/(p)
        y <- (1-u)^(2/3)
        #         y <- ((p-q)/p - u)^(2/3)
        l <- (ell[q+1] + ell[q+2]) / 2
        if (q==0) {
                L <- Inf 
        } else {
                L <- (ell[q] + ell[q+1]) / 2
        }
        #         L <- ell[q+1]
        constrainedFit <- nls(x ~ m + alpha*y, algorithm="port",
                              start=c(m=x[1], alpha=(x[1]-x[2])/(y[1]-y[2])),
#                              start=c(m=x[1], alpha=-6),
                              lower=c(m=l,alpha=-Inf),
                              upper=c(m=L,alpha=Inf)
        )
        return(mean(resid(constrainedFit)^2))
}

set.seed(265)
p <- 50
n <- 50
N <- floor(sqrt(p))
lam <- c(10,4)
k <- length(lam)
max_k <- 10
#DF <- data.frame(iter=rep(1:iter, each=max_k+1), q=rep(0:max_k+1, iter))
iter <- 100
k_hat <- rep(0,iter)
for (i in 1:iter) {
#         Sigma <- diag(c(lam, rep(0,p-k))+rep(1,p))
        Sigma <- diag(c(lam, rep(0,p-k)) + runif(p,0.5,1.5))
        S <- rWishart(1,n,Sigma)
        S <- S[,,1]/n
        ell <- eigen(S,TRUE,TRUE)$values
        mse <- rep(0, max_k+1)
        for (q in 0:max_k) {
                mse[q+1] <- advancedSqrtFitMSE(q, ell, N)
        }
        #plot(0:max_k, log(mse))
        k_hat[i] <- which.min(mse)-1
}
hist(k_hat)
