optimization <- function(x, psi, d){
  psi_pa1 <- vector()
  
  coef <- -lm(log(psi[t,])~., data.frame(x^2, x))$coefficients
  a <- coef[2:(1+d)]
  b <- coef[(2+d):length(coef)]
  
  psi_pa1 <- c(b/(-2*a), 1/(2*a))
  
  return(psi_pa1)
}
