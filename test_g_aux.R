
d <- 2
psi_pa <- matrix(runif(4*d), nrow = 2, byrow = TRUE)
t <- 1
ini_cov <- diag(runif(d), d, d)
A <- diag(runif(d), d, d)
B <- diag(runif(d), d, d)
C <- diag(runif(d), d, d)
D <- diag(runif(d), d, d)
x <- runif(d)
y <- runif(d)
n = 2


test_g_aux <- function(y, x, t, psi_pa, n){  
  expected_output <- dmvn(y, C%*%x, D, log = TRUE) + dmvn(as.vector(A%*%x), psi_pa[t, 1:d], B + diag(psi_pa[t, (d+1):(d+d)], d, d), log = TRUE)-
    dmvn(x, psi_pa[t, 1:d], diag(psi_pa[t, (d+1):(d+d)], d, d), log = TRUE)
  
  result <- g(y, x) + psi_tilda(x, psi_pa, t, n) - psi_t(x, psi_pa, t, n)

tolerance <- 1e-6
is_close_enough <- abs(result - expected_output) < tolerance
if (!is_close_enough) {
  stop("Test failed: Incorrect output from optimization function.")
}

return(0)
  
}
