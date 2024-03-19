iAPF1 <- function(time_step, L, N, Time, d) {
  #initialization for time t = 1...L to obtain psi1, w1, X1
  L <- L
  output <- Init(time_step, L)
  X <- output[[1]] 
  w <- output[[2]]
  psi <- output[[3]]
  
  #choose psi from 1 to 3L/4
  psi_final <- psi[complete.cases(psi), ][1:ceiling(3*L/4),]
  
  #Iteration
  #obtain 
  if(L < Time){
    for(time_step in seq(2*L,Time,L)){
      
      #I didn't include any resampling in this step
      #Run iAPF with the initial distribution we defined 
      output2 <- run_iAPF(time_step, w, X, L, N, Time, d)
      
      X <- output2[[1]]
      w <- output2[[2]]
      psi <- output2[[3]]
      gap_matrix <- matrix(NA, nrow = L/2, ncol = 2*d)
      
      #take psi from kL+L/4+1 to kL+3L/4+1
      psi_final <- rbind(psi_final, gap_matrix, psi[complete.cases(psi), ][ceiling(L/4 + 1):ceiling(3*L/4),])
    }  
  }
  return(psi_final)
}

iAPF2 <- function(psi_final, time_step, L1, L) {
  #t = 1, . . . , L/2 to obtain X2, w2 and psi2
  psi_final <- psi_final
  L1 <- L1
  L <- L
  
  output <- Init(time_step, L1) #pass
  X <- output[[1]] 
  w <- output[[2]]
  psi <- output[[3]]    
  
  ####Algorithm####
  
  #start from L/2 to 3L/2 to obtain psi2
  count = 0
  for(time_step in seq(ceiling(3*L/2),Time, L)){
    #I didn't include any resampling in this step
    #Run iAPF with the initial distribution we defined 
    output2 <- run_iAPF(time_step, w, X, L, N, Time, d)
    
    #smoothing particles
    X <- output2[[1]]
    w <- output2[[2]]
    psi <- output2[[3]]
    
    if(time_step < Time){
      psi_final[ceiling(3*L/4+1+count*L):ceiling(3*L/4+count*L+ L/2),] <- psi[complete.cases(psi), ][ceiling(L/4 + 1):ceiling(3*L/4),]
    }
    
    if(time_step == Time){
      psi_final <- rbind(psi_final, psi[ceiling(nrow(psi_final) + 1):Time,])
    }
    count = count + 1
  }
  
  if(time_step != Time){
    time_step <- Time
    
    output2 <- run_iAPF(time_step, w, X, L1, N, Time, d)
    
    #smoothing particles
    X <- output2[[1]]
    w <- output2[[2]]
    psi <- output2[[3]]
    psi_final <- rbind(psi_final, psi[complete.cases(psi), ][ceiling(L1/2 + 1):L1,])
  }
  
  return(psi_final)
}

quasi_online_main <- function(Lag, N, Time){
  psi_final1 <- iAPF1(Lag, Lag)
  psi_final <- iAPF2(psi_final1, ceiling(Lag/2), ceiling(Lag/2), Lag)
  
  output2 <- smoothing_APF(psi_final, N, Time)
  X <- output2[[1]]
  w <- output2[[2]]
  Z <- output2[[3]]
  
  return(list(X, w, Z))
}
