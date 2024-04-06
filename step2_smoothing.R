smoothing_APF <- function(psi_pa, N, Time,){ 
  #l >= 2
  X_apf <- array(NA, dim = c(Time, N[1], d))
  w_apf <- matrix(NA, Time, N[1])
  Z_apf <- 0
  
  re = 0
  
  #when n = kL, we use the new mu.tilda to initialize particles
  
  X_apf[1,,] <- mu_aux(psi_pa, 1, N, 1)
  for(i in 1:Num){
    w_apf[1,i] <- g_aux_smoo(obs[1], X_apf[1,i,], 1, psi_pa, Time)
  }
  
  for(t in 2:Time){
    
    if(ESS(w_apf[t-1,], t, is.log = TRUE) <= kappa*Num){
      re = re + 1
      output <- residual_resampling(w_apf[t-1,], Num, Z_apf[l])
      mix <- output[[1]]
      Z_apf <- output[[2]]
      
      X_apf <- X_apf[, mix,, drop = FALSE]
      
      for(i in 1:Num){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1, i,], psi_pa, t)
        w_apf[t,i] <- g_aux_smoo(obs[t,], X_apf[t,i,], t, psi_pa, Time) 
        
      }
    }else{
      
      for(i in 1:Num){
        #filtering particles
        X_apf[t,i,] <- f_aux(X_apf[t-1,i,], psi_pa, t) 
        w_apf[t,i] <- w_apf[t-1,i] + g_aux_smoo(obs[t,], X_apf[t,i,], t, psi_pa, Time)
      }
    }
    
  }
  
  mx <- max(w_apf[Time, 1:Num])
  Z_apf <- Z_apf + log(mean(exp(w_apf[Time, 1:Num]-mx))) + mx
  
  return(list(X_apf, w_apf, Z_apf, re))
}
