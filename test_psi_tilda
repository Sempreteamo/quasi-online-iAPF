d <- 2
psi_pa <- matrix(runif(4*d), nrow = 2, byrow = TRUE)
t <- 1
A <- diag(runif(d), d, d)
B <- diag(runif(d), d, d)
xt <- runif(d)

test_psi_tilda <- function(psi_pa, t, d) {
  # Define test data
  # Compute the expected output
  expected_output <- dmvn(as.vector(A%*%xt), psi_pa[t, 1:d], B + diag(psi_pa[t, (d+1):(d+d)], d, d), log = TRUE)
  
  # Test the optimization function
  result <- (-d/2)*log(2*pi) -(1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)]+ diag(B), nrow=d, ncol=d))) +
    (-1/2)*t(A%*%xt - psi_pa[t, 1:d])%*%diag((psi_pa[t, (d+1):(d+d)] + diag(B))^(-1), nrow=d,ncol=d)%*%
    (A%*%xt - psi_pa[t, 1:d]) 
  
  # Check if the result matches the expected output
  tolerance <- 1e-6
  is_close_enough <- abs(result - expected_output) < tolerance
  if (!is_close_enough) {
    stop("Test failed: Incorrect output from optimization function.")
  }
  
  return(0)
}
