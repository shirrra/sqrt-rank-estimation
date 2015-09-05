sqrtFitMSE <- function(q, ell, d=10) {
        p <- length(ell)
        m <- min(ell[ell>0]) # smallest non-zero eigenvalue
        range <- (ell[q+1] - m) / d #we want to consider a range which is 1/d from the total range of eigenvalues
        j <- which.min(ell > (ell[q+1]-range)) #we want to take eigenvalues that are between ell[q+1] and ell[q+1]+range 
        # we should take ell from (q+1) to j
        if (j < q+10) {j <- q+10} #we want to use at least 10 eigenvalues
        x <- ell[(q+1):j]
        u <- ((p-q):(p-j+1))/(p)
        y <- (1-u)^(2/3)
        l <- (ell[q+1] + ell[q+2]) / 2
        if (q==0) {
                L <- Inf 
        } else {
                L <- (ell[q] + ell[q+1]) / 2
        }
        constrainedFit <- lm(y ~ x)
        lm <- -constrainedFit$coefficients[["(Intercept)"]] / constrainedFit$coefficients[["x"]]
        if (lm < l || lm > L) {
                b <- sum(y)/(sum(l-x))
                return(sum((y-b*x)^2)/length(x))
        }
        return(mean(resid(constrainedFit)^2))        
}

kEst <- function(ell, k_max=10) {
        mse <- rep(0,k_max+1)
        for (q in 0:max_k) {
                mse[q+1] <- sqrtFitMSE(q, ell)
        }
        k <- which.min(mse)-1
}

set.seed(26580)
P <- 2^(5:10)
lam <- c(10,4)
k <- length(lam)
max_k <- 10
iter <- 100
frac <- rep(0,length(P))
for (i_P in 1:length(P)) {
        p <- P[i_P]
        n <- p
        count <- 0
        for (i_iter in 1:iter) {
                Sigma <- diag(c(lam, rep(0,p-k)) + runif(p,0.5,1.5))
                S <- rWishart(1,n,Sigma)
                S <- S[,,1]/n
                ell <- eigen(S,TRUE,TRUE)$values
                k_hat <- kEst(ell)
                if (k_hat == k) {count <- count+1}
        }
        frac[i_P] <- count/iter
}
plot(log(P,2), frac, ylim=c(0,1), pch=16, main="Prob. of Estimating k Correctly as a Function of p", xlab="P", ylab="Pr(k_hat=k)")
abline(1, 0, col="red", lwd=3)
