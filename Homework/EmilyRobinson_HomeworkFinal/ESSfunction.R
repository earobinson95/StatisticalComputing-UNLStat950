ESS <- function(chain, stop = 0.1, burnin = 0.5, alpha = 0.05){
  if(!is.matrix(chain)) chain = matrix(chain, ncol = 1)
  if(burnin){
    if (burnin < 1){
      burnin = burnin*dim(chain)[1]
    }
  }
  p = dim(chain)[2]
  results = matrix(NA, p, 7)
  colnames(results) = c("mean", "lowerHPD", "upperHPD","se", "sd", "L", "ESS")
  rownames(results) = colnames(chain)
  for(i in 1:p){
    h = chain[-(1:burnin),i]
    L = length(h)
    hbar = mean(h)
    hdev = h - hbar
    hvar = crossprod(hdev)/L
    tau = 1
    k = 1
    repeat{
      rho = crossprod(hdev[-(1:k)],hdev[-((L+1-k):L)])/((L-k)*hvar)
      tau = tau + 2*rho
      if(abs(rho) < stop || k > 1000) break
      k = k + 1
    }
    ESS = L/tau
    h_sort = sort(h)
    I = matrix(NA, nrow = ((L-1)-floor((1-alpha)*(L-1))), ncol = 2)
    for(j in 1: ((L-1)-floor((1-alpha)*(L-1)))){
      I[j,1] = h_sort[j]
      I[j,2] = h_sort[(j + floor((1-alpha)*(L-1)))]
    }
    jstar = which.min(I[,2] - I[,1])
    HPD = I[jstar,]
    results[i,] = c(mean = hbar, lowerHPD = HPD[1], upperHPD = HPD[2], se = sqrt(hvar/ESS), sd = sqrt(hvar), L = L, ESS = ESS)
  }
  return(results)
}

ESS(furSealChain)
