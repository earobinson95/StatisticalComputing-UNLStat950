baseball <- read.csv("data/baseball.dat", sep="")
summary(baseball)
head(baseball)

# ---------------------------------------------------------------------------------
# a -------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

calcAIC <- function(candidate, X, y){
  n   = dim(X)[1]
  Xc  = as.matrix(cbind(rep(1,n), X[,candidate]))
  pc  = dim(Xc)[2]
  fit = lm.fit(Xc,y)
  return(n*log(crossprod(fit$residuals)/n)+2*(pc+1))
}

localSearch <- function(objectiveFn, candidate, nSteps = 20, reStarts = 5, minimize = FALSE, ...){
  p = length(candidate)
  objectiveValues = rep(NA, (nSteps+1)*reStarts)
  finalCandidates = matrix(0,reStarts,p)
  finalObjective  = rep(NA, reStarts)
  it = 0
  
  for(s in 1:reStarts){
    candidate = (runif(p)<0.5)
    currentObjective = objectiveFn(candidate,...)
    it = it+1
    objectiveValues[it] = NA
    it = it+1
    objectiveValues[it] = currentObjective
    
    for(i in 1:nSteps){
      bestCandidate = candidate
      bestObjective = currentObjective
      for(j in 1:p){
        for(k in j:p){
          newCandidate = candidate
          newCandidate[j] = !newCandidate[j]
          if(k == j){newCandidate[k] = newCandidate[k]
        }else{newCandidate[k] = !newCandidate[k]}
          newObjective = objectiveFn(newCandidate,...)
          
          if(minimize){
            if(newObjective < bestObjective){
              bestObjective = newObjective
              bestCandidate = newCandidate
            }
          } else{ #maximize
            if(newObjective > bestObjective){
              bestObjective = newObjective
              bestCandidate = newCandidate
            }
          }
        }
      }
      currentObjective = bestObjective
      candidate = bestCandidate
      it = it+1
      objectiveValues[it] = currentObjective
    }
    finalCandidates[s,candidate] = 1
    finalObjective[s] = bestObjective
  }
  return(list(it = it,
              finalObjective = finalObjective,
              finalCandidates = finalCandidates,
              objectiveValues = objectiveValues))
}

X <- baseball[,-1]
y <- log(baseball$salary)
candidate <- rep(1,dim(X)[2])
results   <- localSearch(calcAIC, candidate, nSteps = 20, reStarts = 5, minimize = TRUE, X = X, y = y)

# ---------------------------------------------------------------------------------
# b -------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

numsteps <- function(R, S, p, k){R*S*choose(p,k)}

# ---------------------------------------------------------------------------------
# c -------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

numsteps(R = 5, S = 20, p = 27, k = 2)

# ---------------------------------------------------------------------------------
# d -------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------

# Smallest AIC for each of the five restarts
minAIC <- results$finalObjective
minAIC

# Overall best AIC
bestAIC <- min(minAIC)
bestAIC

# ---------------------------------------------------------------------------------
# e -------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
results_data = data.frame(cum_step = seq(1,results$it,1), negAIC = -results$objectiveValues)
plot(results_data$cum_step, results_data$negAIC, type = "l", 
     xlab = "Cumulative Steps", ylab = "Negative AIC", ylim = c(360, 420))
