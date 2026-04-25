priors<- function(x, pi_1=0.2, pi_2=0.5, pi_3 = 1, first_own_lag_prior_mean=NULL, lambda_pr_mean=NULL, lambda_pr_covmat=NULL, Jeffrey=TRUE){
  
  priors <- list()
  
  setup <- x$setup
  yt <- x$data
  m <- setup$m
  p <- setup$p
  d <- setup$d
  xt <- setup$xt
  
  Sigma_AR <- diag(0,m)
  
  for (i in 1:m){
    
    y <- yt[,i]
    N = length(y)-p
    
    Y <- y[-c(1:p)]
    W <- embed(y, dimension = p+1)[, -1]
    X <- xt[-c(1:p), ,drop=F]
    Q <- embed(xt, dimension = p+1)[, -(1:d)]
    
    Z <- cbind(W,X)
    Gamma_hat = solve(crossprod(Z,Z),crossprod(Z,Y))
    U = Y-Z%*%Gamma_hat
    sigma2 <- crossprod(U,U)/(nrow(Z)-ncol(Z))
    Sigma_AR[i,i] <- sigma2
  }
  
  V <- lapply(1:p, function(x) matrix(0, m, m))
  sigma <- sqrt(diag(Sigma_AR))
  
  for (l in 1:p) {
    for (i in 1:m) {
      for (j in 1:m) {
        if (i == j) {
          V[[l]][i,j] <- (pi_1/(l^pi_3))^2
        } else {
          V[[l]][i,j] <- ((pi_1*pi_2*sigma[i])/(l^pi_3*sigma[j]))^2
        }
      }
    }
  }
  V_mat <- do.call(cbind, V)
  gamma_d_pr_covmat <- diag(c(t(V_mat)))
  
  if (is.null(first_own_lag_prior_mean)) first_own_lag_prior_mean <- rep(0,m)
  
  mat <- matrix(0, nrow = m*p, ncol = m)
  for (i in 1:m){
    mat[i,i] <- first_own_lag_prior_mean[i]
  }
  gamma_d_pr_mean = c(mat)
  if(isFALSE(Jeffrey)){
    v_=m+2
    S_ = (m_0-m-1)*setup$Psi_OLS
    priors$S_ <- S_
    priors$v_ <- v_
  }
  
  priors$gamma_d_pr_mean <- gamma_d_pr_mean
  priors$gamma_d_pr_covmat <- gamma_d_pr_covmat
  
  priors$lambda_pr_mean <- lambda_pr_mean
  priors$lambda_pr_covmat <- lambda_pr_covmat
  
  priors$Jeffrey <- Jeffrey
  priors$Sigma_AR <- Sigma_AR
  x$priors <- priors
  
  return(x)
}
