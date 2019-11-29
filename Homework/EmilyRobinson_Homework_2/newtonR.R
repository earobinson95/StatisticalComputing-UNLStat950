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
  
  return(list(x = xNext, se = sqrt(diag(solve(-evalFinal$der2))), value = evalFinal$value, convergence = (Conv < relConvCrit), results = results[1:t,]))
}

ans = newtonR(logLikeSurvival, c(0.5,0.25,5), X = model.matrix(~leukemia$group), y = leukemia$remissiontime, w = as.numeric(!leukemia$censored), relConvCrit = 1.e-14)
results <- cbind(ans$x, ans$se)
greeks = c(alpha='\u03b1', tau='\u03c4', sigma='\u03c3', sigmaSq='\u03c3\u00B2', beta0='\u03b2\u2080', beta1 = '\u03b2\u2081', gamma='\u03b3')
colnames(results) = c("Estimate", "StdErr")
rownames(results) = c(greeks['alpha'], greeks['beta0'], greeks['beta1'])
results
