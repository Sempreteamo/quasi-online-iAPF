#ess resampling

#x = w[t-1,]

ESS <- function(x, t, is.log=FALSE){
  if(is.log) {
    mx <- max(x)
    s <- sum(exp(x - mx))
    ess <- 1/sum((exp(x - mx)/s)^2)
  }else{
    s <- sum(x)
    ess <- 1/sum((x/s)^2) 
  }
  return(ess)  
}
