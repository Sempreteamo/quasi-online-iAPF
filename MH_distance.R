#x = w[Time,]

MH_distance <- function(x, Time, X){
  mx <- max(x) 
  w_ <- exp(x-mx)/sum(exp(x - mx))
  
  #Mahalanobis distance
  for(specific in 1:Time){
    weighted_mean <- colSums(w_*X[specific,,])
    output_s <- Smoothing(specific)
    fks_mean <- output_s[[1]]
    fks_cov <- output_s[[2]]
    MH_distance[specific] <- mahalanobis(weighted_mean, fks_mean, fks_cov)/d
  }
  
  plot(MH_distance)
}
