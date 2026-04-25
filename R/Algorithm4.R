Algorithm4 <- function(x, iter, warmup, H, x_pred, Jeffrey=FALSE){
  #### Algorithm 4 (steady-state BVAR)
  setup <- x$setup
  priors <- x$priors
  Y <- setup$Y
  X <- setup$X
  W <- setup$W
  Q <- setup$Q
  N <- setup$N
  m  <- setup$m
  p  <- setup$p
  d  <- setup$d
  
  Gamma_d_OLS <- setup$Gamma_d_OLS
  Lambda_OLS <- setup$Lambda_OLS
  
  gamma_d_lbar <-  priors$gamma_d_pr_mean
  Sigma_d_lbar <- priors$gamma_d_pr_covmat
  lambda_lbar <- priors$lambda_pr_mean
  Sigma_lambda_lbar <- priors$lambda_pr_covmat
  
  if (isFALSE(Jeffrey)){
    v_ = priors$v_
    S_ = priors$S_
  }
  
  Psi <- vector(mode = "list", length = iter+1)
  gamma_d <- vector(mode = "list", length = iter+1)
  lambda <- vector(mode = "list", length = iter+1)
  
  R <- iter-warmup
  Burnin <- warmup
  forecasts_array <- array(NA, dim = c(H, m, R))
  
  ############### START ###############
  
  #### SELECT STARTING VALUES (OLS ESTIMATES) ####
  
  gamma_d[[1]] <- c(Gamma_d_OLS) #gamma_d^{0}
  lambda[[1]]  <- c(Lambda_OLS) #lambda_d^{0}
  
  for (j in 2:(Burnin+R+1)){
    
    Lambda = matrix(lambda[[j-1]],m,d)
    Gamma_d = matrix(gamma_d[[j-1]],m*p,m)
    
    ############   1   ################
    ############ EQ 29 ################
    U = Y - X%*%t(Lambda) - (W-Q%*%(diag(p) %x% t(Lambda))) %*% Gamma_d
    S = crossprod(U)
    N = nrow(U)
    if (isFALSE(Jeffrey)){
      Psi[[j]] = LaplacesDemon::rinvwishart(N+v_, S+S_)
    } else {
      Psi[[j]] = LaplacesDemon::rinvwishart(N, S)  
    }
    ############   2   ################
    ############ EQ 30 ################
    Y_Lambda = Y-X%*%t(Lambda)
    W_Lambda = (W-Q%*%(diag(p)%x%t(Lambda)))
    
    Sigma_d_bar = solve((solve(Sigma_d_lbar)+solve(Psi[[j]]) %x% (crossprod(W_Lambda))))
    
    gamma_d_bar = Sigma_d_bar %*% (solve(Sigma_d_lbar)%*%gamma_d_lbar + c(t(W_Lambda)%*%Y_Lambda%*%solve(Psi[[j]])))
    
    gamma_d[[j]] = MASS::mvrnorm(1, gamma_d_bar, Sigma_d_bar)
    
    ############   3   ################
    ############ EQ 31 ################
    Gamma_d <- matrix(gamma_d[[j]],m*p,m)
    A_list <- vector("list", p)
    for (i in 1:p) {
      rows_idx <- ((i - 1) * m + 1):(i * m)
      A_list[[i]] <- t(Gamma_d[rows_idx, ])
    }
    blocks <- list(diag(m * d))
    for (i in 1:p) {
      blocks[[i + 1]] <- diag(d) %x% t(A_list[[i]])
    }
    F_prime <- do.call(cbind, blocks)
    F <- t(F_prime)
    
    B = cbind(X,-Q)
    Y_gamma = Y-W%*%Gamma_d
    
    Sigma_lambda_bar = solve(solve(Sigma_lambda_lbar)+t(F) %*% ((crossprod(B))%x%solve(Psi[[j]])) %*% F)
    lambda_bar = Sigma_lambda_bar %*% (solve(Sigma_lambda_lbar)%*%lambda_lbar+t(F)%*%c(solve(Psi[[j]])%*%t(Y_gamma)%*%B))
    
    lambda[[j]] <- MASS::mvrnorm(1, lambda_bar, Sigma_lambda_bar)
    
    ############   4   ################
    if (j > (Burnin+1)) {
    idx <- j - (Burnin + 1)
    Lambda_j <- matrix(lambda[[j]], m, d)
    Gamma_d_j <- matrix(gamma_d[[j]], m*p, m)
    Psi_j <- Psi[[j]]                                    
    
    A_j <- vector("list", p)
    for (i in 1:p) {
      rows_idx <- ((i - 1) * m + 1):(i * m)
      A_j[[i]] <- t(Gamma_d_j[rows_idx,])  # m x m
    }
    
    Y_pred_mat <- matrix(NA, nrow = H, ncol = m)
    
    for (h in 1:H) {
      
      u_t <- MASS::mvrnorm(1, rep(0, m), Psi_j)
      ytilde_t <- x_pred[h, ] %*% t(Lambda_j)
      
      if (h > 1) {
        for (i in 1:min(h - 1, p)) {
          term <- (Y_pred_mat[h - i,] - x_pred[h - i,] %*% t(Lambda_j)) %*% t(A_j[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      
      if (h <= p) {
        for (i in h:p) {
          term <- (Y[N + h - i,] - X[N + h - i,] %*% t(Lambda_j)) %*% t(A_j[[i]])
          ytilde_t <- ytilde_t + term
        }
      }
      
      Y_pred_mat[h, ] <- ytilde_t + u_t
    }
    
    forecasts_array[, ,idx] <- Y_pred_mat
    }
  }
  
  #### The first B (Burnin) draws are discarded as burn-in ###
  keep_idx <- seq(Burnin + 2, Burnin+R+1)
  n_keep <- length(seq(Burnin + 2, Burnin+R+1))
  gamma_d_keep  <- gamma_d[keep_idx]
  lambda_keep <- lambda[keep_idx]
  Psi_keep    <- Psi[keep_idx]
  
  
  ############### FINISH ###############
  
  Gamma_d_draws <- array(unlist(gamma_d_keep), dim = c(m * p, m, n_keep))
  Lambda_draws <- array(unlist(lambda_keep), dim = c(m, d, n_keep))
  Psi_draws <- simplify2array(Psi_keep)
  
  Gamma_d_posterior_mean <- apply(Gamma_d_draws, c(1, 2), mean)
  Lambda_posterior_mean <- apply(Lambda_draws, c(1, 2), mean)
  Psi_posterior_mean <- apply(Psi_draws, c(1, 2), mean)

  res <- list(Gamma_d_draws = Gamma_d_draws,
              Lambda_draws = Lambda_draws,
              Psi_draws = Psi_draws,
              fcst_draws = forecasts_array,
              Gamma_d_posterior_mean = Gamma_d_posterior_mean,
              Lambda_posterior_mean = Lambda_posterior_mean,
              Psi_posterior_mean = Psi_posterior_mean)
  return(res)
}
