library(MASS)

randMat <- function(n, Sigma) {
        p <- dim(Sigma)[1]
        if (p <= n) {
                S <- rWishart(1,n,Sigma)
                S <- S[,,1]/n
        } else {
                X <- mvrnorm(n, rep(0,p), Sigma, )
                S <- 1/n*t(X)%*%X
        }
        return(S)
}

# Fhat <- function(x, l) {
#         p <- length(l)
#         return(sum(l[1:p]<=x) / p)
# }

sqrtMSE1 <- function(k,v,y) {
        b <- v[1]
        x <- (b-v)^1.5
        fit <- lm(y~x)
        return(mean(fit$residuals^2))
}

sqrtMSE2 <- function(k,x,u) {
        y <- (1-u)^(2/3)
        fit <- lm(y~x)
        return(mean(fit$residuals^2))
}

#y = 1-c(m-x)^1.5
#(1-y)^(2/3)=c^(2/3)(m-x)
#y' = 1.5c(m-x)^0.5

KEst1 <- function(l,R) {
        p <- length(l)
        max_k <- 20
        mse <- rep(0,max_k+1)
        for (q in 0:max_k) {
                v <- l[(q+1):(q+R)]
                #y <- sapply(v,Fhat,l=l[(q+1):p])
                u <- ((p-q):(p-q-R+1))/(p-q)
                #y <- sapply(v,Fhat,v)
                mse[q+1] <- sqrtMSE1(q,v,u)
        }
        #plot(0:max_k,log(mse), pch=16)
        return(which.min(mse)-1)
}

KEst2 <- function(l,R) {
        p <- length(l)
        max_k <- 20
        mse <- rep(0,max_k+1)
        for (q in 0:max_k) {
                v <- l[(q+1):(q+R)]
                y <- ((p-q):(p-q-R+1))/(p-q)
                mse[q+1] <- sqrtMSE2(q,v,y)
        }
        #plot(0:max_k,log(mse), pch=16)
        return(which.min(mse[2:(max_k+1)]-mse[1:max_k]))
}

KEst3 <- function(l,R) {
        p <- length(l)
        max_k <- 20
        mse <- rep(0,max_k+1)
        for (q in 0:max_k) {
                v <- l[(q+1):(q+R)]
                #y <- sapply(v,Fhat,l=l[(q+1):p])
                #y <- sapply(v,Fhat,v)
                y <- ((p-q):(p-q-R+1))/(p-q)
                mse[q+1] <- sqrtMSE1(q,v,y)
        }
        #plot(0:max_k,log(mse), pch=16)
#         plot(1:max_k,mse[2:(max_k+1)]-mse[1:max_k], pch=16)
        diff <- mse[2:(max_k+1)]-mse[1:max_k]
        diff_ratio <- diff[1:(max_k-1)]/abs(diff[2:max_k])
        #plot(1:(max_k-1),diff_ratio)
        return(which.min(diff_ratio))
}
