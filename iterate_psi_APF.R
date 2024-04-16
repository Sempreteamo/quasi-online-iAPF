iterate_psi_APF <- function(n, X_apf, Z_apf, w, X, L, N){
  l = 1
  
  while(TRUE){
    
    output <- list()
    
    if(l != 1){
      #print(l)
      #generate filtering particles X_apf for psi the next iteration
      #APF outputs filtering X_apf for the next psi, and smoothing X_apf_s
      #for the final calculation
      
      output <- psi_APF(n, w, X, psi_pa, N[l], L)
      X_apf <- output[[1]]
      w_apf <- output[[2]]
      Z_apf[l] <- output[[3]]
    }
    
    #to speed up the algorithm, I just fix the number of iterations to be k.
    #Here k = 5
    
    if(l <= k ){
      print(l)
      #receive filtering particles X_apf for psi
      psi_pa <- Psi(n, X_apf, N[l], L) 
      
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
