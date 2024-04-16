#N=N[l], x = X_apf[t,,]

Psi <- function(n, x, N, L){
  psi <- matrix(NA, nrow = Time, ncol = N)
  psi_pa1 <- matrix(NA, nrow = Time, ncol = 2*d)
  
  #calculate psi
  for(t in n:(n-L+1)){
    if(t == n){
      for(i in 1:N){
        psi[t,i] <- (1 / ((2 * pi)^(d / 2))) * 
          exp(-0.5 * t(x[t,i,] - obs[t,]) %*% (x[t,i,] - obs[t,]))
      }
      
      
    }else{
      for(i in 1:N){
        
        psi[t,i] <- exp(g(obs[t,],x[t,i,]))*(2*pi)^(-d/2)*prod(psi_pa1[t+1, (d+1):(d+d)]+1)^(-1/2)*
          exp(-(1/2)*t(A%*%x[t,i,] - psi_pa1[t+1, 1:d])%*%diag((psi_pa1[t+1, (d+1):(d+d)]+diag(B))^(-1), nrow=d,ncol=d)%*%
                (A%*%x[t,i,] - psi_pa1[t+1, 1:d]) )
        
      }
    }
    
    psi_pa1[t,] <- optimization(x[t,,], psi[t,])
    
    
    #print(psi_pa1[t, 1:d])
    #print(obs[t,])
    
  }
  return(psi_pa1)
  
}
