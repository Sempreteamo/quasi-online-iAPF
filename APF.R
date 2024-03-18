APF <- function(Time, d, n, w, X, psi_pa, l, Z_apf, N, L){ #purely filtering particles
  #l >= 2
  X_apf <- array(NA, dim = c(Time, N[l], d))
  w_apf <- matrix(NA, Time, N[l])
  Z_apf[l] <- 0
  
  #when n = kL, we use the new mu.tilda to initialize particles
  #X_apf[1:(n-L+1),1:N[l],-] <- mu_aux(psi_pa, l, N, n-L+1)
  #for(i in 1:N[l]){
  # w_apf[1:(n-L+1),i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,],n-L+1, psi_pa, n, L) 
  #}
  if(n == L){
    X_apf[n-L+1,1:N[l],] <- mu_aux(psi_pa, l, N, n-L+1)
    for (i in 1:(n-L)) {
      X_apf[i,,] <- X_apf[n-L+1,,]
    }
    for(i in 1:N[l]){
      w_apf[1:(n-L+1),i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,],n-L+1, psi_pa, n, L) 
    }
  }else{
    output <- change_mupsi(n, w, X, psi_pa, n-L+1, N, l, L)
    X_apf[n-L+1,1:N[l],] <- output[[1]]
    for (i in 1:(n-L)) {
      X_apf[i,,] <- X_apf[n-L+1,,]
    }
    sum_ <- output[[2]]
    
    for (i in 1:N[l]){
      #X_apf[1:(n-L+1),i,] <- f_aux(X[n-L, i,], psi_pa, n-L+1)
      #w_apf[1:(n-L+1), i] <- g_aux(obs[n-L+1,], X_apf[n-L+1,i,], n-L+1, psi_pa, n, L)
      
      w_apf[1:(n-L+1), i] <- g_transition(obs[n-L+1,], X_apf[n-L+1,i,],n-L+1, psi_pa, n) + log(sum_)
    }
    
  }
  
  for(t in (n-L+2):n){
    print(t)
    if(ESS(t,w_apf, is.log = TRUE) <= kappa*N[l]){
      
     output <- multi_resampling(w_apf, t, Z_apf, l)
     mix <- output[[1]]
     Z_apf <- output[[2]]
      #mix <- residual(t, w_apf)
      
      for(i in 1:N[l]){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1, mix[i],], psi_pa, t)
        w_apf[t,i] <- g_aux(obs[t,], X_apf[t,i,], t, psi_pa, n, L) 
      }
    }else{
      
      for(i in 1:N[l]){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1,i,], psi_pa, t) 
        w_apf[t,i] <- w_apf[t-1,i] + g_aux(obs[t,], X_apf[t,i,], t, psi_pa, n, L)
      }
    }
    
  }
  mx <- max(w_apf[n, 1:N[l]])
  Z_apf[l] <- Z_apf[l] + log(mean(exp(w_apf[n, 1:N[l]]-mx))) + mx
  
  return(list(X_apf, w_apf, Z_apf))
}
