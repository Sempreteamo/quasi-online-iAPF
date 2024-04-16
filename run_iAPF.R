run_iAPF <- function(n, w, X, L, N){
  output <- init_APF(n, w, X, L, N)
  X_apf <- output[[1]]
  w_apf <- output[[2]]
  Z_apf <- output[[3]]
  
  output2 <- iterate_psi_APF(n, X_apf, Z_apf, w, X, L, N)
  
  #smoothing particles
  X <- output2[[1]]
  w <- output2[[2]]
  psi <- output2[[3]]
  Z_apf <- output2[[4]]
  
  return(list(X, w, psi))
}
