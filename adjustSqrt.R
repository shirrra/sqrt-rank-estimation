# in this progrem we simulate samples from a rank k factor model 
# with noise, and try to recover the rank k.
# if the noise eigenvalues come from a distribution F_p,n, which converges
# to a distribution F_p/n when p,n->infty with p/n->r, then the sample noise eigenvalues distribution
# converges to T_r(F), a generalized MP distribution.
# let [L,R] be the support of T_r(F), then for some z we have that in
# [R-z,R] T_r(F) behaves like c\sqrt(R-x)
# we use this fact to find R in the finite model

# p - problem dimension
# n - number of samples
# r - p/n
# k - signal space dimension
# x - signals strength
# Sigma - population covariance matrix = diag(x,..,x,0,...0), has rank k
# S - sample covariance matrix

# noise eigvals dist: support [0,1], p(x)=6x(1-x)...
# 3x^2-2x^3 = l/(p-1), l\in[0,1,...,p]...

set.seed(265)
iter <- 1000

p <- 150
n <- 250
r <- p/n
R <- (1+sqrt(r))^2
k <- 3
lam <- 1*R
s <- floor(sqrt(p))
elem <- s
#noise <- rep(1,p-k)
parts <- c(0.1,0.8)
vals <- c(0.5,0.6,1)
num <- c(floor(parts[1]*(p-k)),floor(parts[2]*(p-k)))
noise <- c(rep(vals[1],num[1]),rep(vals[2],num[2]),rep(1,p-k-num[1]-num[2]))
Sigma <- diag(c(rep(lam,k),noise))
S <- rWishart(1,n,Sigma)
S <- S[,,1]/n
l <- eigen(S,TRUE,TRUE)$values
#hist(l,20)      

Fhat <- function(x, l, k) {
        dim <- length(l)
        return(sum(l[1:(dim-k)]<=x) / (dim-k))
}
x1 <- seq(0,4,0.01)

max_k <- 10
for (k_est in 0:max_k) {
        a <- Fhat(l[k_est+1],l,k_est)
        b <- l[k_est+1]
        x2<- l[(k_est+1):(k_est+elem)]
        z <- (b-x2)^1.5
        y <- sapply(x2,Fhat,l=l,k=k_est)
        #alpha <- (a-y)%*%z/(z%*%z)
        #plot(x2,sapply(x2,Fhat,l=l,k=k_est),col="blue")
        #lines(x2,a-alpha*(b-x2)^1.5,col="red")
        #mse[k_est+1] <- (a-alpha*z-y)%*%(a-alpha*z-y)
        fit <- lm(y~z)
        mse[k_est+1] <- mean(fit$residuals^2)
}

par(mfrow=c(2,2))
hist(noise,s)
hist(l,s)
plot(1:s,l[1:s])
plot(0:max_k,mse,pch=16)

#parts <- c(0.1,0.8)
#vals <- c(0.5,0.6,1)
#num <- c(floor(parts[1]*(p-k)),floor(parts[2]*(p-k)))
#noise <- c(rep(vals[1],num[1]),rep(vals[2],num[2]),rep(1,p-k-num[1]-num[2]))
#noise <- (seq(0.5,1,length.out=p-k))^0.1

