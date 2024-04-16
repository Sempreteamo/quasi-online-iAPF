#w = w[t-1,] or w[t-1,]
#Z_apf = Z_apf[l]
#N = N[l]
multi_resampling <- function(w, N, Z_apf){
  mx <- max(w)
  w_ <- exp(w - mx)/sum(exp(w - mx))
  Z_apf = Z_apf + log(mean(exp(w - mx))) + mx
  mix <- sample(1:N, N, replace = TRUE, prob = w_)
  return(list(mix, Z_apf))
}


#w = w[t-1,]
#Z_apf = Z_apf[l]
residual_resampling <- function(w, Z_apf){
  mx <- max(w)
  w_ <- exp(w - mx)/sum(exp(w - mx))
  
  Z_apf = Z_apf + log(mean(exp(w - mx))) + mx
  
  Ntm <- as.integer(Num*w_)
  
  mix <- unlist(lapply(1:Num, function(i) {rep(i, Ntm[i])}))
  mr <- Num - sum(Ntm)
  
  w_hat <- w_ - Ntm/Num
  w_hat <- w_hat*Num/mr
  
  mix <- c(sample(1:Num, mr, replace = TRUE, prob = w_hat), mix)
  
  return(list(mix, Z_apf))
  
} 
