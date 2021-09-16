# set up the tables ------------------------------------------------------------------

grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
grid_log_lambda_long <- as.vector(countirt:::grid_log_lambda)
grid_logZ_long <- as.vector(countirt:::grid_log_Z)
grid_cmp_var_long <- as.vector(countirt:::grid_cmp_var)

# log_lambda_from_grid --------------------------------------------------------------
log_lambda_from_grid <- function(mu, nu){
  
  # for numerical stability, set mu > 500 to 500
  # Warning: only works if in comparison to mean, observation is negligably small
  # and density will be 0 anyways; dcmp is not a stable implementation of density
  # if mean > 500 and density is actually relevant and different from 0
  mu <- ifelse(mu > 200, 200, mu)
  
  if (length(nu) == 1) {
    nu_filled <- rep(nu, length(mu))
  } else {
    nu_filled <- nu
  }
  
  log_lambda <- interp_from_grid_v(grid_mus, grid_nus, grid_log_lambda_long,
                                   mu, nu_filled)
  
  return(log_lambda)
}

# lambda_from_grid -----------------------------------------------------------------
lambda_from_grid <- function(mu, nu){
  
  # for numerical stability, set mu > 500 to 500
  # Warning: only works if in comparison to mean, observation is negligably small
  # and density will be 0 anyways; dcmp is not a stable implementation of density
  # if mean > 500 and density is actually relevant and different from 0
  mu <- ifelse(mu > 200, 200, mu)
  
  lambda_grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  lambda_grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  
  grid_log_lambda_long <- as.vector(countirt:::grid_log_lambda)
  
  if (length(nu) == 1) {
    nu_filled <- rep(nu, length(mu))
  } else {
    nu_filled <- nu
  }
  
  log_lambda <- interp_from_grid_v(lambda_grid_mus, lambda_grid_nus, grid_log_lambda_long,
                                   mu, nu_filled)
  
  return(exp(log_lambda))
}

# logZ_from_grid ------------------------------------------------------------------
logZ_from_grid <- function(mu, nu){
  
  mu <- ifelse(mu > 200, 200, mu)
  
  grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  
  grid_logZ_long <- as.vector(grid_log_Z)
  
  if (length(nu) == 1) {
    nu_filled <- rep(nu, length(mu))
  } else {
    nu_filled <- nu
  }
  
  logZ <- interp_from_grid_v(grid_mus, grid_nus, grid_logZ_long, mu, nu_filled)
  
  return(logZ)
}

# get_var_cmp ------------------------------------------------------------------
get_var_cmp <- function(mu, nu){
  # outputs bicubic interpolant of pre-computed variance
  # mu and nu can be vectors of the same length, or
  # if nu is a scalar, then replicate it to the same length as mu
  if(length(nu)==1) nu = rep(nu, length(mu))
  
  # current fineness of the gridded pre-computed values
  mu <- ifelse(mu > 200, 200, mu)
  
  grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  
  grid_cmp_var_long <- as.vector(grid_cmp_var)
  
  # perform bicubic interpolant on logLambda and logZ
  var <- interp_from_grid_v(grid_mus, grid_nus, grid_cmp_var_long, mu, nu)
  return(var)
}

# dcmp --------------------------------------------------------------------------
dcmp <- function(data, mu, nu, logprob = FALSE) {
  # only works for length(data) = length(mu) = length(nu)
  
  # for numerical stability, set mu > 500 to 500
  # Warning: only works if in comparison to mean, observation is negligably small
  # and density will be 0 anyways; dcmp is not a stable implementation of density
  # if mean > 500 and density is actually relevant and different from 0
  mu <- ifelse(mu > 200, 200, mu)
  
  grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  grid_log_lambda_long <- as.vector(countirt:::grid_log_lambda)
  grid_logZ_long <- as.vector(grid_log_Z)
  
  out <- dcmp_cpp(data, mu, nu, logprob,
                  grid_mus, grid_nus, grid_mus, grid_nus,
                  grid_log_lambda_long, grid_logZ_long)
  
  return(out)
}

# truncated quadrature rule --------------------------------------------------------
quad_rule <- function(n_nodes, thres = Inf, prob = 0) {
  weights_and_nodes <- gaussHermiteData(n_nodes)
  # rescale because we are approximating integral for normal density
  weights_and_nodes$x <- weights_and_nodes$x * sqrt(2)
  weights_and_nodes$w <- weights_and_nodes$w / sqrt(pi)
  trunc <- abs(weights_and_nodes$x) < thres & 
    weights_and_nodes$w > prob
  weights_and_nodes$x <- weights_and_nodes$x[trunc]
  weights_and_nodes$w <- weights_and_nodes$w[trunc] 
  return(weights_and_nodes)
}


