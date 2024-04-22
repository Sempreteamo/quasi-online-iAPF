psi_APF <- function(n, w, X, psi_pa, N, L){ #purely filtering particles
  #l >= 2
  X_apf <- array(NA, dim = c(Time, N, d))
  w_apf <- matrix(NA, Time, N)
  Z_apf <- 0
  
  #when n = kL, we use the new mu.tilda to initialize particles
  #X_apf[1:(n-L+1),1:N[l],-] <- mu_aux(psi_pa, l, N, n-L+1)
  #for(i in 1:N[l]){
  # w_apf[1:(n-L+1),i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,],n-L+1, psi_pa, n, L) 
  #}
  if(n == L){
    X_apf[n-L+1,,] <- mu_aux(psi_pa, N, n-L+1)
    for (i in 1:(n-L)) {
      X_apf[i,,] <- X_apf[n-L+1,,]
    }
    for(i in 1:N){
      w_apf[1:(n-L+1),i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,], n-L+1, psi_pa, n, L) 
    }
  }else{
    output <- change_mupsi(X[n-L,,], w[n-L,], psi_pa, n-L+1, N, l)
    X_apf[n-L+1,,] <- output[[1]]
    for (i in 1:(n-L)) {
      X_apf[i,,] <- X_apf[n-L+1,,]
    }
    sum_ <- output[[2]]
    
    for (i in 1:N){
      #X_apf[1:(n-L+1),i,] <- f_aux(X[n-L, i,], psi_pa, n-L+1)
      #w_apf[1:(n-L+1), i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,], n-L+1, psi_pa, n, L)
      
      w_apf[1:(n-L+1), i] <- g_transition(obs[n-L+1,], X_apf[n-L+1,i,],  n-L+1, psi_pa, n) + log(sum_)
    }
    
  }
  
  for(t in (n-L+2):n){
    #print(t)
    if(ESS(w_apf[t-1,], is.log = TRUE) <= kappa*N){
      
      #output <- multi_resampling(w_apf[t-1,], N, Z_apf)
      output <- residual_resampling(w_apf[t-1,], Z_apf)
      mix <- output[[1]]
      Z_apf <- output[[2]]
      
      for(i in 1:N){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1, mix[i],], psi_pa, t)
        w_apf[t,i] <- g_aux(obs[t,], X_apf[t,i,], t, psi_pa, n, L) 
      }
    }else{
      
      for(i in 1:N){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1,i,], psi_pa, t) 
        w_apf[t,i] <- w_apf[t-1,i] + g_aux(obs[t,], X_apf[t,i,], t, psi_pa, n, L)
      }
    }
    
  }
  mx <- max(w_apf[n, ])
  Z_apf <- Z_apf + log(mean(exp(w_apf[n,]-mx))) + mx
  
  return(list(X_apf, w_apf, Z_apf))
}
