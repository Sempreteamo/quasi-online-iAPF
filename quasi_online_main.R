quasi_online_main <- function(Lag, N, Time){
  psi_final1 <- parallel_iAPF(psi_final = 0, Lag, Lag, Lag, is.parallel1 = TRUE)
  psi_final <- parallel_iAPF(psi_final1, ceiling(Lag/2), ceiling(Lag/2), Lag, is.parallel1 = FALSE)
  
  output2 <- smoothing_APF(psi_final, N, Time)
  X <- output2[[1]]
  w <- output2[[2]]
  Z <- output2[[3]]
  
  return(list(X, w, Z))
}
