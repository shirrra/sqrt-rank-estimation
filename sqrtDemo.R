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
# c - p/n
# k - rank
# x - signals strength
# Sigma - population covariance matrix = diag(x,..,x,0,...0), has rank k
# S - sample covariance matrix
# l - eigenvalues of S

# noise eigvals dist: support [0,1], p(x)=6x(1-x)...
# 3x^2-2x^3 = l/(p-1), l\in[0,1,...,p]...

set.seed(265)
iter <- 20
P <- c(64, 128, 256, 512)
#P <- c(512)
c <- 4 #p/n
lam <- c(200,50,10,5)

k <- length(lam)
sqrtKEst <- rep(0,iter)
prob <- rep(0,length(P))

#par(mfrow=c(2,3))
for (i in 1:length(P)) {
        p <- P[i]
        n <- floor(p/c)
        R <- floor(sqrt(p))
        for (j in 1:iter) {
                noise <- runif(p,0.5,1.5)        
                noise <- rep(1,p)
                Sigma <- diag(c(lam, rep(0,p-k))+noise)
                S <- randMat(n,Sigma)
                l <- eigen(S,TRUE,TRUE)$values
                sqrtKEst[j] <- KEst3(l,R)        
        }
        prob[i] <- sum(sqrtKEst==k)/iter
}

plot(log(P),prob)


# par(mfrow=c(2,2))
# hist(noise, s, main="histogram of noise population eigenvalues")
# hist(l, s, main="histogram of empirical
#      eigenvalues")
# abline(v=m, col="red")
# plot(1:s, l[1:s], main="sqrt(p) largest eigenvalues")
# plot(0:max_k,log(mse), pch=16, 
#      main="mse for sqrt function match for different values of k")
# 
# print(sum(k_hat==k)/iter)
