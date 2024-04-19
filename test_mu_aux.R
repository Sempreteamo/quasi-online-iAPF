d <- 2
psi_pa <- matrix(runif(4*d), nrow = 2, byrow = TRUE)
N <- 100
t <- 1
ini_cov <- diag(runif(d), d, d)
x1 <- runif(d)

test_mu_aux <- function(psi_pa, N, t, d) {
  # Define test data
  # Compute the expected output
  expected_output <- dmvn(x1, rep(0, d), ini_cov)*dmvn(x1, psi_pa[t, 1:d], diag(psi_pa[t, (d+1):(d+d)], d, d))/
    dmvn(rep(0, d), psi_pa[t, 1:d], ini_cov + diag(psi_pa[t, (d+1):(d+d)], d, d))
  
  # Test the optimization function
  result <- dmvn(x1, (psi_pa[t, (d+1):(d+d)]^(-1) + diag(ini_cov)^(-1))^(-1)*
                   (psi_pa[t, (d+1):(d+d)]^(-1)*psi_pa[t,1:d]), 
                 diag((psi_pa[t, (d+1):(d+d)]^(-1) + diag(ini_cov)^(-1))^(-1), nrow=d, ncol = d))
  # Check if the result matches the expected output
  tolerance <- 1e-6
  is_close_enough <- abs(result - expected_output) < tolerance
  if (!is_close_enough) {
    stop("Test failed: Incorrect output from optimization function.")
  }
  
  return(0)
}
