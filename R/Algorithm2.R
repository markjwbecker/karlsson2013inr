Algorithm2 <- function(x, iter, warmup, H, x_pred, Jeffrey=FALSE){
  #### Algorithm 2 (normal-diffuse and independent normal-Wishart)
  setup <- x$setup
  priors <- x$priors
  Y <- setup$Y
  Z <- setup$Z
  N <- setup$N
  m  <- setup$m
  p  <- setup$p
  d  <- setup$d
  
  Gamma_OLS <- setup$Gamma_OLS
  
  gamma_lbar <-  priors$gamma_pr_mean
  Sigma_gamma_lbar <- priors$gamma_pr_covmat
  
  if (isFALSE(Jeffrey)){
    v_ = priors$v_
    S_ = priors$S_
  }
  
  Psi <- vector(mode = "list", length = iter+1)
  gamma <- vector(mode = "list", length = iter+1)
  
  R <- iter-warmup
  Burnin <- warmup
  forecasts_array <- array(NA, dim = c(H, m, R))
  
  ############### START ###############
  
  #### SELECT STARTING VALUE (OLS ESTIMATE) ####
  
  gamma[[1]] <- c(Gamma_OLS) #gamma^{0}
  
  for (j in 2:(Burnin+R+1)){
    
    Gamma = matrix(gamma[[j-1]],(m*p)+(d),m)
    
    ############   1   ################
    ############ EQ 25 ################
    S_bar = crossprod(Y-Z%*%Gamma)
    N = nrow(Y)
    
    if (isFALSE(Jeffrey)){ #independent normal-Wishart
      Psi[[j]] = LaplacesDemon::rinvwishart(N+v_, S_ + S_bar)
    } else { #normal-diffuse
      Psi[[j]] = LaplacesDemon::rinvwishart(N, S_bar)  
    }
    ############   2   ################
    ############ EQ 24 ################
    
    Sigma_gamma_bar = solve(solve(Sigma_gamma_lbar) + solve(Psi[[j]]) %x% crossprod(Z))
    
    gamma_bar = Sigma_gamma_bar %*% (solve(Sigma_gamma_lbar)%*%gamma_lbar + c(t(Z)%*%Y%*%solve(Psi[[j]])))
    
    gamma[[j]] = MASS::mvrnorm(1, gamma_bar, Sigma_gamma_bar)
    
    ############   3   ################
    
    if (j > (Burnin+1)) {
      idx <- j - (Burnin + 1)
      Gamma_j <- matrix(gamma[[j]], (m*p)+(d), m)
      C_j <- t(matrix(t(Gamma_j)[,((m*p)+(1)):((m*p)+(d))],m,d))
      Psi_j <- Psi[[j]]                                    
      
      Phi_j <- vector("list", p)
      for (i in 1:p) {
        rows_idx <- ((i - 1) * m + 1):(i * m)
        Phi_j[[i]] <- t(Gamma_j[rows_idx,])  # m x m
      }
      
      Y_pred_mat <- matrix(NA, nrow = H, ncol = m)
      
      for (h in 1:H) {
        
        u_t <- MASS::mvrnorm(1, rep(0, m), Psi_j)
        ytilde_t <- x_pred[h, ] %*% C_j
        
        if (h > 1) {
          for (i in 1:min(h - 1, p)) {
            term <- Y_pred_mat[h - i,] %*% t(Phi_j[[i]]) #Phi' = A
            ytilde_t <- ytilde_t + term
          }
        }
        
        if (h <= p) {
          for (i in h:p) {
            term <- Y[N + h - i,] %*% t(Phi_j[[i]]) #Phi' = A
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
  gamma_keep  <- gamma[keep_idx]
  Psi_keep    <- Psi[keep_idx]
  
  ############### FINISH ###############
  
  Gamma_draws <- array(unlist(gamma_keep), dim = c((m*p)+(d), m, n_keep))
  Psi_draws <- simplify2array(Psi_keep)
  
  Gamma_posterior_mean <- apply(Gamma_draws, c(1, 2), mean)
  Psi_posterior_mean <- apply(Psi_draws, c(1, 2), mean)
  
  res <- list(Gamma_draws = Gamma_draws,
              Psi_draws = Psi_draws,
              fcst_draws = forecasts_array,
              Gamma_posterior_mean = Gamma_posterior_mean,
              Psi_posterior_mean = Psi_posterior_mean)
  return(res)
}
