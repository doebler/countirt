
# log_lambda_from_grid --------------------------------------------------------------
log_lambda_from_grid <- function(mu, nu){
  
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
  
  log_lambda <- log_lambda_from_grid(mu, nu)
  
  return(exp(log_lambda))
}

# logZ_from_grid ------------------------------------------------------------------
logZ_from_grid <- function(mu, nu){
  
  mu <- ifelse(mu > 200, 200, mu)
  
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
  
  # perform bicubic interpolant on logLambda and logZ
  var <- interp_from_grid_v(grid_mus, grid_nus, grid_cmp_var_long, mu, nu)
  return(var)
}

# dcmp --------------------------------------------------------------------------
dcmp <- function(data, mu, nu, logprob = FALSE) {
  # only works for length(data) = length(mu) = length(nu)
  # TODO besser vektorisieren und die truncation in c++ machen und hier
  # mit min_mu und max_mu arbeiten wie auch in meinen gradienten
  
  mu <- ifelse(mu > 200, 200, mu)
  
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


