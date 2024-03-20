optimization <- function(x, psi, d){
  coef <- -lm(log(psi[t,])~., data.frame(x^2, x))$coefficients
  a <- coef[2:(1+d)]
  b <- coef[(2+d):length(coef)]
  
  psi_pa1[t,] <- c(solve(-2 * diag(a)) %*% b, diag(solve(2 * diag(a))))
  
  return(psi_pa1)
}
