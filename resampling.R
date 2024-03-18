multi_resampling <- function(w_apf, t, Z_apf, l){
  mx <- max(w_apf[t-1,])
  w_ <- exp(w_apf[t-1,]-mx)/sum(exp(w_apf[t-1, ] - mx))
  Z_apf[l] = Z_apf[l] + log(mean(exp(w_apf[t-1,]-mx))) + mx
  mix <- sample(1:N[l],N[l], replace = TRUE, prob = w_)
  return(list(mix, Z_apf))
}

residual <- function(t, w, Num, Z_apf, l){
  mx <- max(w[t-1,])
  w_ <- exp(w[t-1,] - mx)/sum(exp(w[t-1,] - mx))
  
  Z_apf[l] = Z_apf[l] + log(mean(exp(w_apf[t-1,]-mx))) + mx
  
  Ntm <- as.integer(Num*w_)
  
  mix <- unlist(lapply(1:Num, function(i) {rep(i, Ntm[i])}))
  mr <- Num - sum(Ntm)
  
  w_hat <- w_ - Ntm/Num
  w_hat <- w_hat*Num/mr
  
  mix <- c(sample(1:Num, mr, replace = TRUE, prob = w_hat), mix)
  
  return(list(mix, Z_apf))
  
} 
