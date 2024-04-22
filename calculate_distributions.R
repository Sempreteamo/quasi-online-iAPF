
Num_apf <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l])))/mean(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l]))))
}

#trans prob
#N(; Ax, B)
f <- function(x){
  return (rnorm(d, 0, diag(B)) + A%*%x)   #trans prob
}

#obs prob  
#N(; Cx, D)
g <- function(y, x){  
  return ((-d/2)*log(2*pi) + (-1/2)*t(y-x)%*%(y-x)) #obs prob  C%*%x = x
}

#twisted mu
mu_aux <- function(psi_pa, N, t){  
  return(rmvn(N, (psi_pa[t, (d+1):(d+d)]^(-1) + diag(ini_cov)^(-1))^(-1)*
                (psi_pa[t, (d+1):(d+d)]^(-1)*psi_pa[t,1:d]), 
              diag((psi_pa[t, (d+1):(d+d)]^(-1) + diag(ini_cov)^(-1))^(-1), nrow=d, ncol = d)))
}

#twisted g
g_aux <- function(y, x, t, psi_pa, n, L){  
  if(t == (n-L+1)){
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) - (d/2)*log(2*pi)-(1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d))) - 
             (1/2)*t(-psi_pa[t, 1:d])%*%
             diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
             (-psi_pa[t, 1:d]) - psi_t(x, psi_pa, t, n))  #initialisation of g = t=1 or t=L?
  }else{
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n))  #g_2:T 
  }
}

#t=n-L+1 in APF for n!=L
g_transition <- function(y, x, t, psi_pa, n){
  return(g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n))
}

g_aux_smoo <- function(y, x, t, psi_pa, n){  
  if(t == 1){
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) -(d/2)*log((2*pi))-(1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d))) - 
             (1/2)*t(-psi_pa[t, 1:d])%*%
             diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
             (-psi_pa[t, 1:d]) - psi_t(x, psi_pa, t, n))  #initialisation of g = t=1 or t=L?
  }else{
    return(g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n))  #g_2:T 
  }
}

#twisted f
f_aux <- function(x, psi_pa, t){
  return(diag(((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1), nrow=d,ncol=d)%*%
           (A%*%x + diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d]) +
           rnorm(d, 0, ((psi_pa[t, (d+1):(d+d)])^(-1)+1)^(-1)))  #f_2:T 
}

#psi.tilda_t = f (xt, ψt+1); psi.tilda_n = 1
psi_tilda <- function(x, psi_pa, t, n){  #from 0 to T. 0,T = 1 
  if(t == n){
    psi_tilda <- 0
  }else{   #psi_pa_t = psi_t
    psi_tilda <- (-1/2)*log(2*pi) -(1/2)*log(det(diag(psi_pa[t+1, (d+1):(d+d)]+ diag(B), nrow=d, ncol=d))) +
      (-1/2)*t(A%*%x - psi_pa[t+1, 1:d])%*%diag((psi_pa[t+1, (d+1):(d+d)] + diag(B))^(-1), nrow=d,ncol=d)%*%
      (A%*%x-psi_pa[t+1, 1:d]) 
  }
  return(psi_tilda)
}

#ψt(xt) = N (xt; m_t, Σ_t), m_t, Σ_t obtained in Psi function
psi_t <- function(x, psi_pa, t, n){ 
  if(t == (n + 1)){
    psi_t <- 0
  }else{
    psi_t <- -(d/2)*log(2*pi) -(1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)], nrow=d,ncol=d))) +						
      (-1/2)*t(x-psi_pa[t, 1:d])%*%diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%						
      (x-psi_pa[t, 1:d])
  }
  return(psi_t)
}

#the initial distribution mu at time kL
# ∑W_kL(?)^i*f(X_{k−1}^i,(k−1)L, ·), i in 1:N
#sample particles X_{kL-L+1} using the particles at time X_(k-1)L

#x = w[n-L,]
#mus = X[n-L,,]
change_mu <- function(x, mus){
  sam <- matrix(0, Num, d)
  mx <- max(x)
  #cat('mx=',mx)
  w_ <- exp(x - mx)/sum(exp(x - mx))
  #cat('w_=',w_)
  s <- sample(1:Num, size = Num, replace = TRUE, prob = w_) 
  #cat('s=',s)
  for(i in 1:Num){
    sam[i,] <- rnorm(d, 0, diag(B)) + A%*%mus[s[i],]
  }
  return(list(sam, w_))
}

#the initial distribution mu.tilda at time kL
#the particles and weights I use are exactly the same as above, and I use
#f.tilda to generate particles instead of f

#x = w[n-L,]
#mus = X[n-L,,]
change_mupsi <- function(mus, x, psi_pa, t, N, l){
  sam <- matrix(0, Num, d)
  w_adj <- vector()
  mx <- max(x)
  w_ <- exp(x - mx)/sum(exp(x - mx)) #normalized weights of w[n-L,]
  #here we need to adjust the weights 
  
  for(i in 1:Num){
    w_adj[i] <- x[i] + (-1/2)*t(A%*%mus[i,] - psi_pa[t,1:d])%*%diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
      (A%*%mus[i,] - psi_pa[t,1:d])
  }
  
  mx <- max(w_adj)
  w_tilda <- exp(w_adj - mx)/sum(exp(w_adj - mx)) #adjusted normalized weights
  s <- sample(1:Num, size = Num, replace = TRUE, prob = w_tilda) 
  
  #use the adjusted normalized weights to generate particles
  for(i in 1:N){
    sam[i,] <- diag((psi_pa[t, (d+1):(d+d)]^(-1)+1)^(-1), nrow=d,ncol=d)%*%
      (diag(psi_pa[t, (d+1):(d+d)]^(-1), nrow=d,ncol=d)%*%psi_pa[t,1:d] + A%*%mus[s[i],]) + 
      rnorm(d, 0, (psi_pa[t, (d+1):(d+d)]^(-1)+1)^(-1))
  }
  
  #calculate the normalizing constant
  sum_ <- 0
  for(i in 1:N){
    sum_ = sum_ + exp(w_[i])*
      exp((-1/2)*t(A%*%mus[i,] - psi_pa[t,1:d])%*%diag((psi_pa[t, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
            (A%*%mus[i,] - psi_pa[t,1:d]))
  }
  #psi_tilda_0 = sum(w_i*z_c^i)
  sum_ <- (2*pi)^(-1/2)*(det(diag(psi_pa[t, (d+1):(d+d)]+1, nrow=d,ncol=d)))^(-1/2)*sum_
  
  return(list(sam, sum_))
}
