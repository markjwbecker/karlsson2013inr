restrict_Gamma_d <- function(x, restriction_matrix) {
  
  m <- x$setup$m
  p <- x$setup$p
  
  if (!all(dim(restriction_matrix) == c(m * p, m))) {
    stop("restriction_matrix must have dimension (m*p x m)")
  }
  x$setup$restriction_matrix <- restriction_matrix
  
  zero_indices <- which(c(restriction_matrix) == 0)
  
  if (!is.null(x$priors$gamma_d_pr_covmat)) {
    diag(x$priors$gamma_d_pr_covmat)[zero_indices] <- 0.0000001
  } else {
    warning("gamma_d_pr_covmat not found in priors")
  }
  return(x)
}
