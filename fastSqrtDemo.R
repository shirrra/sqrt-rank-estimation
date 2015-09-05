set.seed(265)
iter <- 9
P <- c(64, 128, 256, 512)
P <- 512
p=512
c <- 0.25 #p/n
lam <- c(200,50)

k <- length(lam)
sqrtKEst <- rep(0,iter)
#prob <- rep(0,length(P))

par(mfrow=c(3,3))
        #p <- P[i]
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
        prob <- sum(sqrtKEst==k)/iter


#plot(log(P),prob)
print(prob)

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
