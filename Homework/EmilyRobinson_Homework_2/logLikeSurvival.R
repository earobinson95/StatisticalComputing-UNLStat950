leukemia <- read.csv("leukemia.dat", sep="")
head(leukemia)

logLikeSurvival <- function(theta, der = 0, X, y, w){
  
  p     <- length(theta)
  alpha <- theta[1]
  beta  <- theta[2:3]
  
  value <- sum(w*log(y^alpha*exp(X%*%beta)) - y^alpha*exp(X%*%beta) + w*log(alpha/y))
  if(der == 0) return(value)
  
  der1      <- matrix(NA, nrow = 3, ncol = 1)  
  der1[1]   <- sum(w*log(y) - y^alpha*log(y)*exp(X%*%beta) + w/alpha)
  for (j in 2:p){
    der1[j]  <- sum(w*X[,j-1] - y^alpha*X[,j-1]*exp(X%*%beta))
  }
  if(der ==1) return(list(value = value, der1 = der1))
  
  der2      <- matrix(NA, nrow = 3, ncol = 3)
  der2[1,1] <- sum(-y^alpha*log(y)^2*exp(X%*%beta)-w/alpha^2)
  for (j in 2:p){
    der2[1,j] <- der2[j,1] <- sum(-y^alpha*X[,j-1]*log(y)*exp(X%*%beta))
    for (k in 2:p){
      der2[j,k] = der2[k,j] = sum(-y^alpha*X[,j-1]*exp(X%*%beta))
    }
  }
  return(list(value = value, der1 = der1, der2 = der2))
  
}

trial <- logLikeSurvival(c(2,2,2), der = 2, X = model.matrix(~leukemia$group), y = leukemia$remissiontime, w = as.numeric(!leukemia$censored))
trial

