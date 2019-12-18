changePointGibbs <- function(x, gamma, nSamples = 10^4){
  x = as.matrix(x)
  N = dim(x)[1]
  # Initial parameters
  theta   = sample(1:N-1,1)
  alpha   = rgamma(1, gamma[2], gamma[3])
  lambda1 = rgamma(1, gamma[1], alpha)
  lambda2 = rgamma(1, gamma[1], alpha)
  
  chain = matrix(NA, nSamples + 1, 4)
  rownames(chain) = 0:nSamples
  
  
  greeks = c(theta = "\u03B8", alpha='\u03b1', lambda1='\u03BB\u2081', lambda2 = '\u03BB\u2082')
  colnames(chain) = c(greeks['theta'], greeks['alpha'], greeks['lambda1'], greeks['lambda2'])
  chain[1,] = c(theta, alpha, lambda1, lambda2)
  
  for(s in 1:nSamples){
    ptheta <- rep(NA, (N-1))
    for(n in 1:(N-1)){
      ptheta[n] = sum(x[1:n])*log(lambda1) - n*lambda1 + sum(x[(n+1):N])*log(lambda2) - (N-n)*lambda2
    }
    prob = exp(ptheta - max(ptheta))/sum(exp(ptheta - max(ptheta)))
    
    theta = which.max(rmultinom(1,1,prob))
    alpha   = rgamma(1, 2*gamma[1] + gamma[2], lambda1 + lambda2 + gamma[3])
    lambda1 = rgamma(1, gamma[1]+sum(x[1:theta]), theta + alpha)
    lambda2 = rgamma(1, gamma[1]+sum(x[(theta+1):N]), N-theta+alpha)
    
    chain[s+1,] = c(theta, alpha, lambda1, lambda2)
  }
  return(chain)
}

coal <- read.table("coal.dat", header = T)
coalChain = changePointGibbs(coal$disasters, gamma = c(3,10,10), nSamples = 10^4)

par(mfrow = c(2,2))
plot(coalChain[,1], xlab = "Samples", ylab = expression(theta))
plot(coalChain[,2], xlab = "Samples", ylab = expression(alpha), type = "l")
plot(coalChain[,3], xlab = "Samples", ylab = expression(lambda[1]), type = "l")
plot(coalChain[,4], xlab = "Samples", ylab = expression(lambda[2]), type = "l")
par(mfrow =c(1,1))

ESS(coalChain)

