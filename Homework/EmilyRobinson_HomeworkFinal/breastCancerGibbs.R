# Import Dataset
breastCancer <- read.table("breastcancer.dat", header = T)
head(breastCancer)
summary(breastCancer)

breastCancerGibbs <- function(X, t, w, alpha, lambda, nSamples = 10^4){
  p = dim(X)[2]
  n = dim(X)[1]
  
  # Initial parameters
  theta = rep(NA, p)
  for(j in 1:p){
    theta[j] = rgamma(1, alpha[j], lambda[j])
  }
  chain = matrix(NA, nSamples + 1, p)
  rownames(chain) = 0:nSamples
  colnames(chain) = c(paste("theta", 1:p, sep = "_"))
  chain[1,] = theta
  for(s in 1:nSamples){
    
    kappaIndex <- matrix(NA, nrow = n, ncol = p)
    for(i in 1:n){
      for(j in 1:p){
        kappaIndex[i,j] <- theta[j]^X[i,j]
      }
    }
    kappa <- apply(kappaIndex, 1, prod)
    
    for(j in 1:p){
      index    = which(X[,j] == 1)
      theta[j] = rgamma(1,sum(w[index]*X[index,j])+alpha[j], lambda[j]+sum(t[index]*kappa[index]/kappaIndex[index,j]))
    }
    chain[s+1,] = theta
  }
  return(chain)
}

startTime = Sys.time()
breastCancerChain = breastCancerGibbs(X = model.matrix(~breastCancer$treatment), 
                                      t = breastCancer$recurtime, 
                                      w = !breastCancer$censored, 
                                      alpha = c(2,2),
                                      lambda = c(60,1),
                                      nSamples = 10^4)
endTime = Sys.time()
endTime - startTime

par(mfrow = c(2,1))
plot(breastCancerChain[,1], xlab = "Samples", ylab = expression(theta[1]), type = "l")
plot(breastCancerChain[,2], xlab = "Samples", ylab = expression(theta[2]), type = "l")
par(mfrow =c(1,1))
ESS(breastCancerChain)

# Parallel Computing
require(parallel)
require(doParallel)

nCores = detectCores()
nChains = 5
ncl = min(nChains, nCores - 1)

registerDoParallel(ncl)
startTime = Sys.time()
chains = foreach(c = 1:nChains) %dopar%{
  breastCancerChain = breastCancerGibbs(X = model.matrix(~breastCancer$treatment), 
                                        t = breastCancer$recurtime, 
                                        w = !breastCancer$censored, 
                                        alpha = c(2,2),
                                        lambda = c(60,1),
                                        nSamples = 10^4)
}
endTime = Sys.time()
endTime - startTime
stopImplicitCluster()

par(mfrow = c(2,1))
plot(chains[[1]][1:1000,1], xlab = "Samples", ylab = expression(theta[1]), col = "black", type = "l")
lines(chains[[2]][1:1000,1], col = "gray28")
lines(chains[[3]][1:1000,1], col = "gray45")
lines(chains[[4]][1:1000,1], col = "gray70")
lines(chains[[5]][1:1000,1], col = "gray87")

plot(breastCancerChain[1:1000,2], xlab = "Samples", ylab = expression(theta[2]), col = "black", type = "l")
lines(chains[[2]][1:1000,2], col = "gray28")
lines(chains[[3]][1:1000,2], col = "gray45")
lines(chains[[4]][1:1000,2], col = "gray70")
lines(chains[[5]][1:1000,2], col = "gray87")
par(mfrow =c(1,1))

head(breastCancerChain)
breastCancerChain2 <- cbind(breastCancerChain, 'kappaC' = 1/breastCancerChain[,1], 'kappaT' = 1/(breastCancerChain[,1]*breastCancerChain[,2]))
round(ESS(breastCancerChain2),3)
