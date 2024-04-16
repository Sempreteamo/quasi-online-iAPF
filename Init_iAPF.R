#Below is the iAPF algorithm
#First we do the initialization
####init_APF####
init_APF <- function(n, w, X, L, N){  #pure filtering
  Z_apf = 0
  X_init <- array(NA, dim = c(Time, Num, d))
  w_init <- matrix(NA, Time, N)
  
  #when n = L, we use the mu as the initialization distribution mu; 
  #when n = kL, we use the new distribution
  #I tried to modify the previous paths from 1 to n-L+1
  #X_init[1:(n-L+1),,] <- rnorm(N*d)  
  #for(i in 1:N){
  # w_init[1:(n-L+1),i] <- g(obs[n-L+1,], X_init[n-L+1,i,])  
  #}
  
  if(n == L){
    X_init[n-L+1,,] <- rnorm(N*d)  
    for (i in 1:(n-L)) {
      X_init[i,,] <- X_init[n-L+1,,]
    }
    for(i in 1:N){
      w_init[1:(n-L+1),i] <- g(obs[n-L+1,], X_init[n-L+1,i,])  
    }
  }else{
    cat('n=',n,'L=',L)
    cat('w[n-L]=', w[n-L,])
    output <- change_mu(w[n-L,], X[n-L,,])
    X_init[n-L+1,,] <- output[[1]]
    for (i in 1:(n-L)) {
      X_init[i,,] <- X_init[n-L+1,,]
    }
    w_ <- output[[2]]
    print('pass')
    
    
    for (i in 1:Num){
      # X_init[1:(n-L+1),i,] <- f(X[n-L,i,])
      w_init[1:(n-L+1), i] <- g(obs[n-L+1,], X_init[n-L+1,i,]) 
    }
    
  }
  
  for(t in (n-L+2):n){
    
    if(ESS(w_init[t-1,], is.log=TRUE) <= kappa*N){
      
      output <- residual_resampling(w_init[t-1,], Z_apf)
      mix <- output[[1]]
      Z_apf <- output[[2]]
      # at the initialization stage, we want filtering particles for psi
      
      for(i in 1:N){
        X_init[t,i,] <- f(X_init[t-1,mix[i],]) 
        w_init[t,i] <- g(obs[t,], X_init[t,i,])  
      }
      
    }else{
      for(i in 1:N){
        X_init[t,i,] <- f(X_init[t-1,i,]) 
        w_init[t,i] <- w_init[t-1,i] + g(obs[t,], X_init[t,i,])  
      }
    }
  }
  
  mx <- max(w_init[n, 1:N])
  Z_apf <- Z_apf + log(mean(exp(w_init[n, 1:N]-mx))) + mx
  return(list(X_init = X_init, w_init = w_init, Z_apf))
}

