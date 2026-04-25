setup <- function(x, ...) UseMethod("setup")

setup.bvar <- function(x, p, deterministic=c("constant", "constant_and_dummy", "constant_and_trend"), dummy=NULL, trend=NULL) {
  
  yt <- x$data
  N = nrow(yt)-p
  m = ncol(yt)
  
  deterministic <- match.arg(deterministic)
  
  if (deterministic == "constant") {
    xt <- cbind(rep(1, nrow(yt)))
    d <- 1
  } else if (deterministic == "constant_and_dummy") {
    xt <- cbind(rep(1, nrow(yt)), dummy)
    d <- 2
  } else {
    trend <- 1:nrow(yt)
    xt <- cbind(rep(1, nrow(yt)), trend)
    d <- 2
  }
  
  Y <- yt[-c(1:p), ]
  W <- embed(yt, dimension = p+1)[, -(1:m)]
  X <- xt[-c(1:p), ,drop=F]
  Q <- embed(xt, dimension = p+1)[, -(1:d), drop=F]
  Z <- cbind(W,X)
  
  Gamma_OLS = solve(crossprod(Z),crossprod(Z,Y))
  U = Y-Z%*%Gamma_OLS
  Psi_OLS <- crossprod(U)/(N-m*p-d)
  
  if (d == 1) {
    C_hat <- Gamma_OLS[(m*p+1):(m*p+d),]
  } else {
    C_hat <- t(Gamma_OLS[(m*p+1):(m*p+d),])
  }
  A <- vector("list", p)
  for (i in 1:p) {
    rows_idx <- ((i - 1) * m + 1):(i * m)
    A[[i]] <- matrix(t(Gamma_OLS[rows_idx, ]), m, m)
  }
  A_L <- diag(m)
  for (i in 1:p) {
    A_L <- A_L - A[[i]]
  }
  Lambda_OLS <- solve(A_L, C_hat)
  Gamma_d_OLS = Gamma_OLS[1:(m*p),]
  x$setup <- list(N=N,
                  m=m,
                  p=p,
                  Y=Y,
                  X=X,
                  W=W,
                  Q=Q,
                  d=d,
                  dummy=dummy,
                  trend=trend,
                  Gamma_OLS=Gamma_OLS,
                  Gamma_d_OLS=Gamma_d_OLS,
                  Lambda_OLS=Lambda_OLS,
                  Psi_OLS=Psi_OLS,
                  xt=xt)
  return(x)
}
