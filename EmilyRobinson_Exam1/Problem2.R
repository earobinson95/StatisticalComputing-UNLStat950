wine <- read.csv("data/wine.dat", sep="")
head(wine)

calcTWGSS <- function(candidate, X){
    K = length(unique(candidate))
    p = dim(X)[2]
    SS <- matrix(NA, K, p)
    
    for(k in 1:K){
      for(j in 1:p){
        x_ij <- X[candidate == k,j]
        xbar_kj <- mean(x_ij)
        SS[k,j] <- sum((x_ij - xbar_kj)^2)
      }
    }
  return(sum(SS))
}

X <- wine[,-1]
y <- wine[,1]
calcTWGSS(y, X)


geneticAlgo <- function(objectiveFn, candidate, G = 99, P = 25, muRate = 0.01, maximum = TRUE, K = 2, ...){
  p = length(candidate)
  multiplier = ifelse(maximum, 1, -1)
  candidates = matrix(sample(1:K, P*p, replace = T),P, p)
  offspring = matrix(NA, P, p)
  objectiveValuesMat = matrix(NA, P, G+1)
  bestObjectiveValue = -Inf
  objectiveValuesMat[,1] = apply(candidates, 1, function(c){objectiveFn(c,...)})
  for(i in 1:P){
    if(multiplier*objectiveValuesMat[i,1] > bestObjectiveValue){
      bestCandidate = (1:p)[candidates[i,]]
      bestObjectiveValue = multiplier*objectiveValuesMat[i,1]
    }
  }
  for(g in 1:G){
    fitness = 2*rank(multiplier*objectiveValuesMat[,g])/(P*(P+1))
    for(i in 1:P){
      parents = rbind(candidates[sample(P,1,prob = fitness),],
                      candidates[sample(P,1),])
      crossover = runif(p)<0.5
      offspring[i,] = parents[1,]
      offspring[i,crossover] = parents[2,crossover]
      mutations = runif(p) < muRate
      offspring[i,mutations] = sample(1:K, length(mutations[mutations == T]), replace = T)
    }
    candidates = offspring
    objectiveValuesMat[,g+1] = apply(candidates,1,function(c){objectiveFn(c,...)})
    for(i in 1:P){
      if(multiplier*objectiveValuesMat[i,g+1] > bestObjectiveValue){
        bestCandidate = (1:p)[candidates[i,]]
        bestObjectiveValue = multiplier*objectiveValuesMat[i,g+1]
      }
    }
  }
  return(list(bestValue = multiplier*bestObjectiveValue, bestCandidate = bestCandidate, objectiveValues = objectiveValuesMat))
}

X <- wine[,-1]
candidate <- rep(1,dim(X)[1])
results   <- geneticAlgo(calcTWGSS, candidate, G = 99, P = 25, muRate = 0.01, maximum = FALSE, K = 4, X = X)
results$bestValue
results$bestCandidate
length(results$bestCandidate)
cbind(wine[,1], results$bestCandidate)
