generate_obs <- function(){
  #set.seed(123)
  X_true[1,] <- rnorm(d, 0, diag(B))     
  for(t in 2:Time){ 
    #set.seed(123)#observations
    X_true[t,] <- rnorm(d, 0, diag(B)) + A%*%X_true[t-1,]  #t(rmvn(d) + A%*%x)
  }
  #set.seed(123)
  return(matrix(rnorm(Time*d, X_true%*%C, D), ncol = d))
}
