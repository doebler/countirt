
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

# marg_ll2 --------------------------------------------------------------------------

marg_ll2 <- function(data, item_params, weights_and_nodes, family, fix_disps = NULL,
                     fix_alphas = NULL, interp_method = "bicubic",
                     same_disps = FALSE, same_alphas = FALSE) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  if (family == "poisson") {
    if (is.null(fix_alphas)) {
      if (same_alphas) {
        alpha <- item_params[grepl("alpha", names(item_params))]
        alphas <- rep(alpha, n_items) 
      } else {
        alphas <- item_params[grepl("alpha", names(item_params))]
      }
    } else {
      alphas <- fix_alphas
    }
  } else if (family == "cmp") {
    if (is.null(fix_alphas)) {
      if (same_alphas) {
        alpha <- item_params[grepl("alpha", names(item_params))]
        alphas <- rep(alpha, n_items)
      } else {
        alphas <- item_params[grepl("alpha", names(item_params))]
      }
    } else {
      alphas <- fix_alphas
    }
    if (is.null(fix_disps)) {
      if (same_disps) {
        log_disp <- item_params[grepl("log_disp", names(item_params))]
        disps <- rep(exp(log_disp), n_items)
      } else {
        log_disps <- item_params[grepl("log_disp", names(item_params))]
        disps <- exp(log_disps)
      }
    } else {
      disps <- fix_disps
    }
  }
  
  
  if (family == "poisson") {
    # function to compute integral with quadrature over
    f <- function(z, data, alphas, deltas) {
      out <- 0
      for (j in 1:n_items) {
        lambda <- exp(alphas[j] * z + deltas[j])
        out <- out + (dpois(data[,j], lambda, log = TRUE))
      }
      return(exp(out))
    }
    
    marg_prob <- numeric(n_persons)
    for (i in 1:n_persons) {
      marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                             data = data[i, , drop = FALSE], 
                             alphas = alphas, deltas = deltas)
    }
    ll <- sum(log(marg_prob))
  } else if (family == "cmp") {
    if (interp_method == "bicubic") {
      ll <- marg_ll_cpp(data = as.matrix(data),
                        alphas = alphas,
                        deltas = deltas, 
                        disps = disps, 
                        nodes = weights_and_nodes$x,
                        weights = weights_and_nodes$w,
                        grid_mus = grid_mus,  
                        grid_nus = grid_nus, 
                        grid_logZ_long = grid_logZ_long,
                        grid_log_lambda_long = grid_log_lambda_long,
                        max_mu = 150,
                        min_mu = 0.001)
    } else {
      # then interpolation method is linear
      # here we don't cap mu, so we extrapolate beyond grid values
      ll <- marg_ll_cpp_lininterp(data = as.matrix(data),
                                  alphas = alphas,
                                  deltas = deltas, 
                                  disps = disps,
                                  nodes = weights_and_nodes$x,
                                  weights = weights_and_nodes$w,
                                  grid_mus = grid_mus,  
                                  grid_nus = grid_nus, 
                                  grid_logZ_long = grid_logZ_long,
                                  grid_log_lambda_long = grid_log_lambda_long)
    }
    
  }
  return(ll)
}

