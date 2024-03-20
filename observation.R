Obs <- function(d, Time,){
  #set.seed(123)
  X_true[1,] <- rnorm(d)     
  for(t in 2:Time){ 
    #set.seed(123)#observations
    X_true[t,] <- rnorm(d) + A%*%X_true[t-1,]  #t(rmvn(d) + A%*%x)
  }
  #set.seed(123)
  return(matrix(rnorm(Time*d, X_true, 1), ncol = d))
}
