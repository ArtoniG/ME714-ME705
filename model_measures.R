model_measures <- function(result)
{
  
  X <- model.matrix(result)
  logLikr <- logLik(result)
  n <- nrow(X)
  p <- ncol(X) #+ 1
  AICe <- AIC(result)
  BICe <- BIC(result)
  AICce <- AICe + 2*p*(p+1)/(n-p-1)
  SABICe <- -2*logLikr + p*log((n+2)/24)
  HQCICe <- -2*logLikr + 2*p*log(log(n))
  m2loglik <- -2*logLikr
  resultE <- c(AICe,BICe,AICce,SABICe,HQCICe,m2loglik)
  names(resultE)<- c("AIC","BIC","AICc","SABIC","HQCIC","-2log.lik")
  return(resultE)
  cat(estatísticac,resultE,"\n")
}
