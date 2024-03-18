
psi_APF <- function(n, X_apf, Z_apf, w, X, L, N, Time, d){
  l = 1
  
  while(TRUE){
    
    output <- list()
    
    if(l != 1){
      print(l)
      #generate filtering particles X_apf for psi the next iteration
      #APF outputs filtering X_apf for the next psi, and smoothing X_apf_s
      #for the final calculation
      
      output <- APF(Time, d, n, w, X, psi_pa, l, Z_apf, N, L)
      X_apf <- output[[1]]
      w_apf <- output[[2]]
      Z_apf <- output[[3]]
    }
    
    #to speed up the algorithm, I just fix the number of iterations to be k.
    #Here k = 5
    
    if(l <= k ){
      
      #receive filtering particles X_apf for psi
      psi_pa <- Psi(l, n, X_apf, N, L, Time) 
      
      if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z_apf[max(l-k,1):l])){  
        N[l+1] <- 2*N[l]
        
      }else{
        N[l+1] <- N[l]
      }
      
      l <- l+1
    }else break
  }
  
  #output psi
  return(list(X_apf, w_apf, psi_pa, Z_apf[l]))
}
