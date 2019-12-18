# Import Dataset
breastCancer <- read.table("breastcancer.dat", header = T)
head(breastCancer)
summary(breastCancer)

# Create Objective Function
logLikeCancer <- function(theta, der = 0, X, t, w){
  p = dim(X)[2]
  n = dim(X)[1]
  
  kappaIndex <- matrix(NA, nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      kappaIndex[i,j] <- theta[j]^X[i,j]
    }
  }
  kappa <- apply(kappaIndex, 1, prod)
  value <- sum(w*log(kappa)-kappa*t)
  if(der == 0) return(value)
  
  der1      <- matrix(NA, nrow = p, ncol = 1)
  for (j in 1:p){
    der1[j]  <- sum((w/kappa - t)*X[,j]*kappa/theta[j])
  }
  if(der ==1) return(list(value = value, der1 = der1))
  
  der2      <- matrix(NA, nrow = p, ncol = p)
  for (j in 1:p){
    for (m in 1:p){
      if(j == m){
        der2[j,m] = sum((-w/kappa^2)*(X[,j]*kappa/theta[j])^2)
      } else {
        der2[j,m] = der2[m,j] = sum((w/kappa-t)*(X[,j]*X[,m]*kappa/(theta[j]*theta[m])) 
                                    + (-w/kappa^2)*(X[,j]*X[,m]*kappa^2)/(theta[j]*theta[m]))
      }
    }
  }
  return(list(value = value, der1 = der1, der2 = der2))
  
}
logLikeCancer(theta = results[,1], der = 2, X = model.matrix(~breastCancer$treatment), t = breastCancer$recurtime, w = !breastCancer$censored)
logLikeCancer(theta = c(0.008,1.2), der = 2, X = model.matrix(~breastCancer$treatment), t = breastCancer$recurtime, w = !breastCancer$censored)

# Check Objective Function
theta1 = 1
theta2 = 1
delta = 0.00000000001

evalTheta       <- logLikeCancer(theta = c(theta1,theta2), der = 2, X = model.matrix(~breastCancer$treatment), t = breastCancer$recurtime, w = !breastCancer$censored)
evalTheta1Delta <- logLikeCancer(theta = c(theta1+delta, theta2), der = 2, X = model.matrix(~breastCancer$treatment), t = breastCancer$recurtime, w = !breastCancer$censored)
evalTheta2Delta <- logLikeCancer(theta = c(theta1, theta2+delta), der = 2, X = model.matrix(~breastCancer$treatment), t = breastCancer$recurtime, w = !breastCancer$censored)

# Check score function
evalTheta$der1[1]/((evalTheta1Delta$value - evalTheta$value)/delta)
evalTheta$der1[2]/((evalTheta2Delta$value - evalTheta$value)/delta)

# Check Hessian
evalTheta$der2[1,1]/((evalTheta1Delta$der1[1] - evalTheta$der1[1])/delta)
evalTheta$der2[1,2]/((evalTheta2Delta$der1[1] - evalTheta$der1[1])/delta)

evalTheta$der2[2,1]/((evalTheta1Delta$der1[2] - evalTheta$der1[2])/delta)
evalTheta$der2[2,2]/((evalTheta2Delta$der1[2] - evalTheta$der1[2])/delta)

# Newton Function
newtonR <- function(f, xInit, maxIt = 20, relConvCrit = 1.e-10, ...){
  p = length(xInit)
  results = matrix(NA, maxIt, p+2)
  colnames(results) = c("value", paste("x", 1:p, sep = ""), "Conv")
  
  xCurrent = xInit
  for(t in 1:maxIt){
    evalF = f(xCurrent, der = 2, ...)
    results[t, "value"] = evalF$value
    results[t, 1+(1:p)] = xCurrent
    xNext = xCurrent - solve(evalF$der2, evalF$der1)
    Conv = sqrt(crossprod(xNext - xCurrent))/(sqrt(crossprod(xCurrent))+relConvCrit)
    results[t, "Conv"] = Conv
    if(Conv < relConvCrit) break
    xCurrent = xNext
  }
  
  evalFinal <- f(xNext, der = 2, ...)
  
  return(list(x = xNext, se = sqrt(diag(-solve(evalFinal$der2))), value = evalFinal$value, convergence = (Conv < relConvCrit), results = results[1:t,]))
}

ans = newtonR(logLikeCancer, c(0.01,1), X = model.matrix(~breastCancer$treatment), t = breastCancer$recurtime, w = !breastCancer$censored, relConvCrit = 1.e-16)
results <- cbind(ans$x, ans$se)
greeks = c(theta1 = "θ1", theta2 = "θ2")
colnames(results) = c("Estimate", "StdErr")
rownames(results) = c(greeks['theta1'], greeks['theta2'])
round(results,3)
