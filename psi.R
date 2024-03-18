Psi <- function(l, n, X_apf, N, L, Time){
  psi <- matrix(NA, nrow = Time, ncol = N[l])
  psi_pa1 <- matrix(NA, nrow = Time, ncol = 2*d)
  
  #calculate psi
  for(t in n:(n-L+1)){
    if(t == n){
      for(i in 1:N[l]){
        psi[t,i] <- (1 / ((2 * pi)^(d / 2))) * 
          exp(-0.5 * t(X_apf[t,i,] - obs[t,]) %*% (X_apf[t,i,] - obs[t,]))
      }
      
      
    }else{
      for(i in 1:N[l]){
        
        psi[t,i] <- exp(g(obs[t,],X_apf[t,i,]))*dmvn(as.vector(A%*%X_apf[t,i,]), 
                                                     psi_pa1[t+1, 1:d], diag(psi_pa1[t+1, (d+1):(d+d)]+1, nrow=d,ncol=d))
      }
    }
    
    psi_pa1 <- optimization(X_apf, t, psi, d)
    
    
    print(psi_pa1[t, ])
    print(obs[t,])
    
  }
  return(psi_pa1)
  
}
