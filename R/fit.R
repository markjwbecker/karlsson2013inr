fit <- function(x, iter = 5000, warmup = 2500, steadystate=FALSE) {
  
  Jeffrey <- x$priors$Jeffrey
    
  if (steadystate) {
    x$fit$SteadyState <- Algorithm4(
      x = x,
      iter = iter,
      warmup = warmup,
      H = x$predict$H,
      x_pred = x$predict$x_pred,
      Jeffrey = Jeffrey
    )
    
    x$posterior_means$Gamma_d <- x$fit$gibbs$Gamma_d_posterior_mean
    x$posterior_means$Lambda <- x$fit$gibbs$Lambda_posterior_mean
    x$posterior_means$Psi <- x$fit$gibbs$Psi_posterior_mean 
    
  } else {
    x$fit$Algorithm2 <- Algorithm2(
      x = x,
      iter = iter,
      warmup = warmup,
      H = x$predict$H,
      x_pred = x$predict$x_pred,
      Jeffrey = Jeffrey
    )
    
    x$posterior_means$Gamma <- x$fit$gibbs$Gamma_posterior_mean
    x$posterior_means$Psi <- x$fit$gibbs$Psi_posterior_mean 
  }
    return(x)
}