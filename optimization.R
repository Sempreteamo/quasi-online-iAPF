optimization <- function(X_apf, t, psi, d){
  data <- cbind(X_apf[t,,]^2, X_apf[t,,])
  coef <- -lm(log(psi[t,])~., as.data.frame(data))$coefficients
  c <- coef[1]
  a <- coef[2:(1+d)]
  b <- coef[(2+d):length(coef)]
  
  psi_pa1[t,] <- c(solve(-2 * diag(a)) %*% b, diag(solve(2 * diag(a))))
  
  return(psi_pa1)
}
