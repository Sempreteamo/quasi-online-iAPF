library(mvtnorm)
library(MASS)
library(FKF)
library(profvis)
library(mvnfast)

##parameters
start = proc.time()
set.seed(72323)
Num <- 500 #total number of particles
N <- vector()
N[1] <- Num
Time = 200
Lag = 8 #lag can be any integers <= Time which is divided by Time
alpha = 0.42
d = 2
k <- 5
tau <- 0.5
kappa = 0.5
A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}
B = C = D = diag(1, nrow = d, ncol = d)
X_true <-  matrix(0, nrow = Time, ncol = d )
obs <-  matrix(0, nrow = Time, ncol = d )
dt <- ct <- matrix(0, d, 1)
Tt <- A
P0 <- Zt <- Ht <- Gt <- diag(1, d, d)
a0 <- rep(0, d)
index <- 1
psi_final <- matrix(NA, nrow = Time, ncol = 2*d)
MH_distance <- vector()
KS_distance <- vector()
X <- array(NA, dim = c(Time, Num, d))
X_ <- array(NA, dim = c(Time, Num, d))
w <- matrix(NA, Time, Num)
Z = 0
Z_apf <- vector()


file_paths <- list.files("C:\\Users\\14775\\Documents\\quasi-online-iapf", full.names = TRUE)

for (file_path in file_paths) {
  source(file_path)
}

obs <- generate_obs()
fkf.obj <- normalizing_constant(a0, P0, dt, ct, Tt, Zt, Ht, Gt, obs)

quasi_online_main <- function(Lag){
  psi_final1 <- parallel_iAPF(psi_final = 0, Lag, Lag, Lag, is.parallel1 = TRUE)
  psi_final <- parallel_iAPF(psi_final1, ceiling(Lag/2), ceiling(Lag/2), Lag, is.parallel1 = FALSE)
  
  output2 <- smoothing_APF(psi_final)
  X <- output2[[1]]
  w <- output2[[2]]
  Z <- output2[[3]]
  
  return(list(X, w, Z))
}


output <- quasi_online_main(Lag)
Z <- output[[3]]

ratio <- log_ratio(Z, fkf.obj)
