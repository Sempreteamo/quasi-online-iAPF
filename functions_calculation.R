#ess resampling
ESS <- function(t,w, is.log=FALSE){
  if(is.log) {
    mx <- max(w[t-1,])
    s <- sum(exp(w[t-1,]-mx))
    ess <- 1/sum((exp(w[t-1,]-mx)/s)^2)
  }else{
    s <- sum(w[t-1,])
    ess <- 1/sum((w[t-1,]/s)^2) 
  }
  return(ess)  
}

residual <- function(t, w){
  mx <- max(w[t-1,])
  w_ <- exp(w[t-1,] - mx)/sum(exp(w[t-1,] - mx))
  
  Ntm <- as.integer(Num*w_)
  
  mix <- unlist(lapply(1:Num, function(i) {rep(i, Ntm[i])}))
  mr <- Num - sum(Ntm)
  
  w_hat <- w_ - Ntm/Num
  w_hat <- w_hat*Num/mr
  
  mix <- c(sample(1:Num, mr, replace = TRUE, prob = w_hat), mix)
  
  return(mix)
  
} 

Num_apf <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l])))/mean(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l]))))
}

#trans prob
#N(; Ax, B)
f <- function(x){
  return (rnorm(d) + A%*%x)   #trans prob
}

#obs prob  
#N(; Cx, D)
g <- function(y, x){  
  return (log((2*pi)^(-d/2)) + (-1/2)*t(y-x)%*%(y-x)) #obs prob  C%*%x = x
}

#twisted mu
mu_aux <- function(psi_pa, l, N, t){  
  return(rmvn(N[l], diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                (diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d]), 
              diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))
}

#twisted g
g_aux <- function(y, x, t, psi_pa, n, L){  
  if(t == (n-L+1)){
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) -(1/2)*log(2*pi)-(1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d))) - 
             (1/2)*t(-psi_pa[t, 1:d])%*%
             diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
             (-psi_pa[t, 1:d]) - psi_t(x, psi_pa, t, n))  #initialisation of g = t=1 or t=L?
  }else{
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n))  #g_2:T 
  }
}

#t=n-L+1 in APF for n!=L
g_transition <- function(y, x, t, psi_pa, n){
  return(g(y, x) + psi_tilda(x, psi_pa, t, n) -
           psi_t(x, psi_pa, t, n))
}
#(1/2)*t(-psi_pa[t, 1:d])%*%
# diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
#(-psi_pa[t, 1:d])

g_aux_smoo <- function(y, x, t, psi_pa, n){  
  if(t == 1){
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) -(1/2)*log((2*pi))-(1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d))) - 
             (1/2)*t(-psi_pa[t, 1:d])%*%
             diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
             (-psi_pa[t, 1:d]) - psi_t(x, psi_pa, t, n))  #initialisation of g = t=1 or t=L?
  }else{
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n))  #g_2:T 
  }
}

#twisted f
f_aux <- function(x, psi_pa, t){
  return(rmvn(1, diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                (A%*%x + diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d]), 
              diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)))  #f_2:T 
}

#psi.tilda_t = f (xt, ψt+1); psi.tilda_n = 1
psi_tilda <- function(x, psi_pa, t, n){  #from 0 to T. 0,T = 1 
  if(t == n){
    psi_tilda <- 0
  }else{   #psi_pa_t = psi_t
    psi_tilda <- (-1/2)*log(2*pi) -(1/2)*log(det(diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d, ncol=d))) +
      (-1/2)*t(A%*%x - psi_pa[t+1, 1:d])%*%diag((psi_pa[t+1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
      (A%*%x-psi_pa[t+1, 1:d]) 
  }
  return(psi_tilda)
}

#ψt(xt) = N (xt; m_t, Σ_t), m_t, Σ_t obtained in Psi function
psi_t <- function(x, psi_pa, t, n){ 
  if(t == (n + 1)){
    psi_t <- 0
  }else{
    psi_t <- -(1/2)*log(2*pi) -(1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)], nrow=d,ncol=d))) +						
      (-1/2)*t(x-psi_pa[t, 1:d])%*%diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%						
      (x-psi_pa[t, 1:d])
  }
  return(psi_t)
}

#the distribution here and the distribution below should be double checked

#the initial distribution mu at time kL
# ∑W_kL(?)^i*f(X_{k−1}^i,(k−1)L, ·), i in 1:N
#I used the weights at time (k-1)L because these are the closest weights, 
#sample particles X_{kL-L+1} using the particles at time X_(k-1)L

change_mu <- function(n, w, X, L){
  sam <- matrix(0,Num,d)
  mx <- max(w[n-L,])
  w_ <- exp(w[n-L,]-mx)/sum(exp(w[n-L,] - mx))
  s <- sample(1:Num, size = Num, replace = TRUE, prob = w_) 
  mus <- X[n-L,,]
  for(i in 1:Num){
    sam[i,] <- rnorm(d) + A%*%mus[s[i],]
  }
  return(list(sam, w_))
}

#the initial distribution mu.tilda at time kL
#the particles and weights I use are exactly the same as above, and I use
#f.tilda to generate particles instead of f
change_mupsi <- function(n, w, X, psi_pa, t, N, l, L){
  sam <- matrix(0,Num,d)
  w_adj <- vector()
  mx <- max(w[n-L,])
  w_ <- exp(w[n-L,]-mx)/sum(exp(w[n-L,] - mx)) #normalized weights of w[n-L,]
  #here we need to adjust the weights 
  
  for(i in 1:Num){
    w_adj[i] <- w[n-L,i] + (-1/2)*t(A%*%X[n-L,i,] - psi_pa[t,1:d])%*%diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
      (A%*%X[n-L,i,] - psi_pa[t,1:d])
  }
  
  mx <- max(w_adj)
  w_tilda <- exp(w_adj-mx)/sum(exp(w_adj - mx)) #adjusted normalized weights
  s <- sample(1:Num, size = Num, replace = TRUE, prob = w_tilda) 
  mus <- X[n-L,,]
  
  #use the adjusted normalized weights to generate particles
  for(i in 1:N[l]){
    sam[i,] <- rmvn(1,  diag((psi_pa[t, (d+1):(d+d)]^(-1)+1)^(-1), nrow=d,ncol=d)%*%
                      (diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d] + A%*%mus[s[i],]),
                    diag((psi_pa[t, (d+1):(d+d)]^(-1)+1)^(-1), nrow=d,ncol=d))
  }
  
  
  #calculate the normalizing constant
  sum_ <- 0
  for(i in 1:N[l]){
    sum_ = sum_ + exp(w_[i])*
      exp((-1/2)*t(A%*%X[n-L,i,] - psi_pa[t,1:d])%*%diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
            (A%*%X[n-L,i,] - psi_pa[t,1:d]))
    
  }
  #psi_tilda_0 = sum(w_i*z_c^i)
  sum_ <- (2*pi)^(-1/2)*(det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d)))^(-1/2)*sum_
  
  return(list(sam, sum_))
}
