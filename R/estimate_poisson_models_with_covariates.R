
# estep_poisson_with_cov --------------------------------------------------------------------

estep_poisson_with_cov <- function(data, item_params, 
                                   p_covariates, i_covariates, 
                                   weights_and_nodes,
                                   i_cov_on = c("alpha", "delta")) {
  
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params)) & 
                          !grepl("beta", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params)) & 
                          !grepl("beta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    # e step for person covariates
    p_covariates <- as.matrix(p_covariates)
    
    PPs <- matrix(
      log(weights_and_nodes$w),
      nrow = nrow(data),
      ncol = length(weights_and_nodes$x),
      byrow = TRUE
    )
    
    for (j in 1:ncol(data)) {
      lambdas <- exp(outer(
        alphas[j] * as.numeric(p_covariates %*% betas_p),
        alphas[j] * weights_and_nodes$x + deltas[j],
        "+"
      ))
      PPs <- PPs + apply(lambdas, 2, function(x){dpois(data[,j], x, log = TRUE)})
    }
  } else if (!is.null(i_covariates)) {
    # e step for item covariates
    
    # distinguish between on which parameters we have the item covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") { # case with item covariates on delta
        i_covariates <- as.matrix(i_covariates)
        
        PPs <- matrix(
          log(weights_and_nodes$w),
          nrow = nrow(data),
          ncol = length(weights_and_nodes$x),
          byrow = TRUE
        )
        
        sum_icov <- as.numeric(i_covariates %*% betas_i)
        for (j in 1:ncol(data)) {
          lambdas <- exp(deltas + alphas[j] * weights_and_nodes$x + sum_icov[j])
          # note that deltas will be just one scalar value in the case of item covariates
          PPs <- PPs + outer(data[,j], lambdas, dpois, log = TRUE)
        }
      } else if (i_cov_on == "alpha") { # case with item covariates on alpha
        i_covariates <- as.matrix(i_covariates)
        
        PPs <- matrix(
          log(weights_and_nodes$w),
          nrow = nrow(data),
          ncol = length(weights_and_nodes$x),
          byrow = TRUE
        )
        
        sum_icov <- as.numeric(i_covariates %*% betas_i)
        for (j in 1:ncol(data)) {
          lambdas <- exp(deltas[j] + alphas * weights_and_nodes$x + weights_and_nodes$x * sum_icov[j])
          # note that deltas will be just one scalar value in the case of item covariates
          PPs <- PPs + outer(data[,j], lambdas, dpois, log = TRUE)
        }
      }
    } else {
      # the alternative is here only that we have covariates on both item parameters
      i_covariates <- as.matrix(i_covariates)
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      
      PPs <- matrix(
        log(weights_and_nodes$w),
        nrow = nrow(data),
        ncol = length(weights_and_nodes$x),
        byrow = TRUE
      )
      
      sum_icov_alpha <- as.numeric(i_covariates %*% betas_i_alpha)
      sum_icov_delta <- as.numeric(i_covariates %*% betas_i_delta)
      for (j in 1:ncol(data)) {
        lambdas <- exp(deltas + alphas * weights_and_nodes$x +
                         weights_and_nodes$x * sum_icov_alpha[j] + 
                         sum_icov_delta[j])
        # note that deltas will be just one scalar value in the case of item covariates
        PPs <- PPs + outer(data[,j], lambdas, dpois, log = TRUE)
      }
    } # end case with covariates on both item parameters
  }
  
  PPs <- exp(PPs)
  
  PPs <- PPs / rowSums(PPs)
  
  # output should be a matrix with N rows and K cols
  return(PPs)
}

# grad_poisson_with_cov -----------------------------------------------------------------------

grad_poisson_with_cov <- function(item_params, PPs, weights_and_nodes, data,
                                  p_covariates, i_covariates,
                                  i_cov_on = c("alpha", "delta")) {
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params)) & 
                          !grepl("beta", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params)) & 
                          !grepl("beta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  N <- nrow(data)
  M <- ncol(data)
  K <- length(weights_and_nodes$x)
  I <- length(betas_i)
  P <- length(betas_p)
  nodes <- weights_and_nodes$x
  
  if (!is.null(p_covariates)) {
    # model with person covariates
    grad_alphas <- numeric(length(alphas))
    grad_deltas <- numeric(length(deltas))
    grad_betas_p <- numeric(length(betas_p))
    p_covariates <- as.matrix(p_covariates)
    for (j in 1:M) {
      for (k in 1:K) {
        lambda <- exp(deltas[j] + alphas[j] * nodes[k] +
                        # alphas[j] * rowSums(t(apply(p_covariates, 1, function(x){x*betas_p}))))
                        as.numeric(p_covariates%*%betas_p * alphas[j]))
        # 2nd line is going to yield a vector of length N which we want so that then
        # our lambda is person specific
        grad_deltas[j] <- grad_deltas[j] + sum((data[,j] - lambda)*PPs[,k])
        grad_alphas[j] <- grad_alphas[j] + sum((nodes[k] +  as.numeric(p_covariates%*%betas_p))*
                                                  #rowSums(t(apply(p_covariates, 1, function(x){x*betas_p}))))*
                                                 (data[,j] - lambda)*PPs[,k])
      }
    }
    for (p in 1:P) {
      for (k in 1:K) {
        for (j in 1:M) {
          lambda <- exp(deltas[j] + alphas[j] * nodes[k] +
                          # alphas[j] * rowSums(t(apply(p_covariates, 1, function(x){x*betas_p}))))
                           as.numeric(p_covariates%*%betas_p * alphas[j]))
          # p_covariates%*%betas_p is going to yield a vector of length N which we want so that then
          # our lambda is person specific
          grad_betas_p[p] <- grad_betas_p[p] + sum(alphas[j]*p_covariates[,p]*(data[,j] - lambda)*PPs[,k])
        }
      }
    }
    out <- c(grad_alphas, grad_deltas, grad_betas_p)
  } else if (!is.null(i_covariates)) {
    # model with item covariates
    
    # distinguish between on which parameters we have the item covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        grad_alphas <- numeric(length(alphas))
        grad_delta <- 0
        grad_betas_i <- numeric(length(betas_i))
        i_covariates <- as.matrix(i_covariates)
        for (j in 1:M) {
          for (k in 1:K) {
            lambda <- exp(deltas + alphas[j] * nodes[k] + sum(as.numeric(betas_i*i_covariates[j,])))
            # note that for item covariates on delta, deltas is a scalar
            grad_alphas[j] <- grad_alphas[j] + sum(nodes[k]*(data[,j] - lambda)*PPs[,k])
            grad_delta <- grad_delta + sum((data[,j] - lambda)*PPs[,k])
          }
        }
        for (c in 1:I) {
          for (k in 1:K) {
            for (j in 1:M) {
              lambda <- exp(deltas + alphas[j] * nodes[k] + sum(as.numeric(betas_i*i_covariates[j,])))
              # note that for item covariates on delta, deltas is a scalar
              grad_betas_i[c] <- grad_betas_i[c] + sum(i_covariates[j,c]*(data[,j] - lambda)*PPs[,k])
            }
          }
        }
        out <- c(grad_alphas, grad_delta, grad_betas_i)
      } else if (i_cov_on == "alpha") {
        grad_alpha <- 0
        grad_deltas <- numeric(length(deltas))
        grad_betas_i <- numeric(length(betas_i))
        i_covariates <- as.matrix(i_covariates)
        for (j in 1:M) {
          for (k in 1:K) {
            lambda <- exp(deltas[j] + alphas * nodes[k] + nodes[k] * sum(as.numeric(betas_i*i_covariates[j,])))
            # note that for item covariates on alpha, alphas is a scalar
            grad_alpha <- grad_alpha + sum(nodes[k]*(data[,j] - lambda)*PPs[,k])
            grad_deltas[j] <- grad_deltas[j] + sum((data[,j] - lambda)*PPs[,k])
          }
        }
        for (c in 1:I) {
          for (k in 1:K) {
            for (j in 1:M) {
              lambda <- exp(deltas[j] + alphas * nodes[k] + nodes[k] * sum(as.numeric(betas_i*i_covariates[j,])))
              # note that for item covariates on alpha, alphas is a scalar
              grad_betas_i[c] <- grad_betas_i[c] + sum(i_covariates[j,c]*nodes[k]*(data[,j] - lambda)*PPs[,k])
            }
          }
        }
        out <- c(grad_alpha, grad_deltas, grad_betas_i)
      }
    } else {
      # the alternative is just that we have item parameters on all parameters
      I <- length(betas_i)/2
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      grad_alpha <- 0
      grad_delta <- 0
      grad_betas_i_alpha <- numeric(length(betas_i)/2)
      grad_betas_i_delta <- numeric(length(betas_i)/2)
      i_covariates <- as.matrix(i_covariates)
      for (j in 1:M) {
        for (k in 1:K) {
          lambda <- exp(deltas + alphas * nodes[k] + 
                          nodes[k] * sum(as.numeric(betas_i_alpha*i_covariates[j,])) +
                          sum(as.numeric(betas_i_delta*i_covariates[j,])))
          # alphas and deltas are covariates when we have covariates on both covariates
          grad_alpha <- grad_alpha + sum(nodes[k]*(data[,j] - lambda)*PPs[,k])
          grad_delta <- grad_delta + sum((data[,j] - lambda)*PPs[,k])
        }
      }
      for (c in 1:I) {
        for (k in 1:K) {
          for (j in 1:M) {
            lambda <- exp(deltas + alphas * nodes[k] + 
                            nodes[k] * sum(as.numeric(betas_i_alpha*i_covariates[j,])) +
                            sum(as.numeric(betas_i_delta*i_covariates[j,])))
            # alphas and deltas are covariates when we have covariates on both covariates
            grad_betas_i_alpha[c] <- grad_betas_i_alpha[c] + 
              sum(i_covariates[j,c]*nodes[k]*(data[,j] - lambda)*PPs[,k])
            grad_betas_i_delta[c] <- grad_betas_i_delta[c] +
              sum(i_covariates[j,c]*(data[,j] - lambda)*PPs[,k])
          }
        }
      }
      out <- c(grad_alpha, grad_delta, grad_betas_i_alpha, grad_betas_i_delta)
    }
    
  } # end case with item covariates

  return(out)
}

# grad_poisson_with_cov_fixalphas --------------------------------------------------------------

grad_poisson_with_cov_fixalphas <- function(item_params, PPs, weights_and_nodes, 
                                   data, p_covariates, i_covariates,
                                   fix_alphas, i_cov_on = "delta") {
  # note: if we fix alphas, then we can only implement item covariates on delta,
  # because we are not estimating alpha at all
  # so we also can't have covariates on both parameters at once
  data <- as.matrix(data)
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  N <- nrow(data)
  M <- ncol(data)
  K <- length(weights_and_nodes$x)
  I <- length(betas_i)
  P <- length(betas_p)
  nodes <- weights_and_nodes$x
  
  if (!is.null(p_covariates)) {
    # model with person covariates
    grad_deltas <- numeric(length(deltas))
    grad_betas_p <- numeric(length(betas_p))
    p_covariates <- as.matrix(p_covariates)
    for (j in 1:M) {
      for (k in 1:K) {
        lambda <- exp(deltas[j] + alphas[j] * nodes[k] + as.numeric(p_covariates%*%betas_p * alphas[j]))
        # p_covariates%*%betas_p is going to yield a vector of length N which we want so that then
        # our lambda is person specific
        grad_deltas[j] <- grad_deltas[j] + sum((data[,j] - lambda)*PPs[,k])
      }
    }
    for (p in 1:P) {
      for (k in 1:K) {
        for (j in 1:M) {
          lambda <- exp(deltas[j] + alphas[j] * nodes[k] + as.numeric(p_covariates%*%betas_p * alphas[j]))
          # p_covariates%*%betas_p is going to yield a vector of length N which we want so that then
          # our lambda is person specific
          grad_betas_p[p] <- grad_betas_i[c] + sum(alphas[j]*p_covariates[,p]*(data[,j] - lambda)*PPs[,k])
        }
      }
    }
    out <- c(grad_deltas, grad_betas_p)
  } else if (!is.null(i_covariates)) {
    # model with item covariates
    grad_delta <- 0
    grad_betas_i <- numeric(length(betas_i))
    i_covariates <- as.matrix(i_covariates)
    for (j in 1:M) {
      for (k in 1:K) {
        lambda <- exp(deltas + alphas[j] * nodes[k] + sum(as.numeric(betas_i*i_covariates[j,])))
        # note that for item covariates, deltas is a scalar
        grad_delta <- grad_delta + sum((data[,j] - lambda)*PPs[,k])
      }
    }
    for (c in 1:I) {
      for (k in 1:K) {
        for (j in 1:M) {
          lambda <- exp(deltas + alphas[j] * nodes[k] + sum(as.numeric(betas_i*i_covariates[j,])))
          # note that for item covariates, deltas is a scalar
          grad_betas_i[c] <- grad_betas_i[c] + sum(i_covariates[j,c]*(data[,j] - lambda)*PPs[,k])
        }
      }
    }
    out <- c(grad_delta, grad_betas_i)
  }
  
  return(out)
}

# grad_poisson_with_cov_samealpha --------------------------------------------------------------

grad_poisson_with_cov_samealpha <- function(item_params, PPs, weights_and_nodes, 
                                            data, p_covariates, i_covariates,
                                            i_cov_on = "delta") {
  # note: same as with fix_alphas, if we want the same alpha across items, it doesn't
  # make sense to include item covariates because they have different values for different
  # items and would necessarily then lead to different alphas
  # so we also can't have covariates on both parameters at once here
  data <- as.matrix(data)
  alphas <- rep(item_params[grepl("alpha", names(item_params))], ncol(data))
  alpha <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  grad_alpha <- 0
  if (is.null(i_covariates)) {
    grad_deltas <- numeric(length(deltas))
  }
  # in the case of item covariates we have just one delta (which is the intercept)
  # in our prediction of item-specific delta_j
  
  N <- nrow(data)
  M <- ncol(data)
  K <- length(weights_and_nodes$x)
  I <- length(betas_i)
  P <- length(betas_p)
  nodes <- weights_and_nodes$x
  
  if (!is.null(p_covariates)) {
    # model with person covariates
    grad_betas_p <- numeric(length(betas_p))
    p_covariates <- as.matrix(p_covariates)
    for (j in 1:M) {
      for (k in 1:K) {
        lambda <- exp(deltas[j] + alphas[j] * nodes[k] + as.numeric(p_covariates%*%betas_p * alphas[j]))
        # p_covariates%*%betas_p is going to yield a vector of length N which we want so that then
        # our lambda is person specific
        grad_deltas[j] <- grad_deltas[j] + sum((data[,j] - lambda)*PPs[,k])
        grad_alpha <- grad_alpha + sum((nodes[k]+p_covariates%*%betas_p)*(data[,j] - lambda)*PPs[,k])
      }
    }
    for (p in 1:P) {
      for (k in 1:K) {
        for (j in 1:M) {
          lambda <- exp(deltas[j] + alphas[j] * nodes[k] + as.numeric(p_covariates%*%betas_p * alphas[j]))
          # p_covariates%*%betas_p is going to yield a vector of length N which we want so that then
          # our lambda is person specific
          grad_betas_p[p] <- grad_betas_i[c] + sum(alphas[j]*p_covariates[,p]*(data[,j] - lambda)*PPs[,k])
        }
      }
    }
    out <- c(grad_alpha, grad_deltas, grad_betas_p)
  } else if (!is.null(i_covariates)) {
    # model with item covariates
    grad_delta <- 0
    grad_betas_i <- numeric(length(betas_i))
    i_covariates <- as.matrix(i_covariates)
    for (j in 1:M) {
      for (k in 1:K) {
        lambda <- exp(deltas + alphas[j] * nodes[k] + sum(as.numeric(betas_i*i_covariates[j,])))
        # note that for item covariates, deltas is a scalar
        grad_alpha <- grad_alpha + sum(nodes[k]*(data[,j] - lambda)*PPs[,k])
        grad_delta <- grad_delta + sum((data[,j] - lambda)*PPs[,k])
      }
    }
    for (c in 1:I) {
      for (k in 1:K) {
        for (j in 1:M) {
          lambda <- exp(deltas + alphas[j] * nodes[k] + sum(as.numeric(betas_i*i_covariates[j,])))
          # note that for item covariates, deltas is a scalar
          grad_betas_i[c] <- grad_betas_i[c] + sum(i_covariates[j,c]*(data[,j] - lambda)*PPs[,k])
        }
      }
    }
    out <- c(grad_alpha, grad_delta, grad_betas_i)
  }
  
  return(out)
}

# ell_poisson_with_cov ----------------------------------------------------------------------------

ell_poisson_with_cov <- function(item_params, PPs, weights_and_nodes, 
                                 data, p_covariates, i_covariates,
                                 i_cov_on = c("alpha", "delta")) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params)) & !grepl("beta", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params)) & !grepl("beta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  N <- nrow(data)
  M <- ncol(data)
  K <- length(weights_and_nodes$x)
  I <- length(betas_i)
  P <- length(betas_p)
  
  nodes <- weights_and_nodes$x
  
  out <- 0
  if (!is.null(p_covariates)) { # we have person covariates
    for (k in 1:K) {
      for (i in 1:N) {
        for(j in 1:M) {
          log_mu <- alphas[j] * nodes[k] + deltas[j]
          for (p in 1:P) {
            log_mu <- log_mu + betas_p[p] * alphas[j] * p_covariates[i,p]
          }
          mu <- exp(log_mu);
          out <- out + dpois(data[i,j], mu, log = TRUE)*PPs[i,k]
        }
      }
    }
  } else if (!is.null(i_covariates)) { # we have item covariates
    
    # then we distinguish between which model parameters are predicted through
    # the item covariates
    
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") { # case of covariates on delta
        for (k in 1:K) {
          for (i in 1:N) {
            for(j in 1:M) {
              log_mu <- deltas + alphas[j] * nodes[k]
              # for item covariates deltas is just one scalar value when we have covariates
              # on delta
              for (c in 1:I) {
                log_mu <- log_mu + betas_i[c] * i_covariates[j,c]
              }
              mu <- exp(log_mu);
              out <- out + dpois(data[i,j], mu, log = TRUE)*PPs[i,k]
            }
          }
        }
      } else if (i_cov_on == "alpha") { # case of covariates on alpha
        for (k in 1:K) {
          for (i in 1:N) {
            for(j in 1:M) {
              log_mu <- deltas[j] + alphas * nodes[k]
              # for item covariates alphas is just one scalar value when we have covariates
              # on alpha
              for (c in 1:I) {
                log_mu <- log_mu + betas_i[c] * i_covariates[j,c] * nodes[k]
              }
              mu <- exp(log_mu);
              out <- out + dpois(data[i,j], mu, log = TRUE)*PPs[i,k]
            }
          }
        }
      }
    } else {
      # for the poisson case, the alternative is just covariates on all item parameters
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      for (k in 1:K) {
        for (i in 1:N) {
          for(j in 1:M) {
            log_mu <- deltas + alphas * nodes[k]
            # for item covariates on both, both alphas and deltas are scalars
            for (c in 1:I) {
              log_mu <- log_mu + betas_i_alphas[c] * i_covariates[j,c] * nodes[k] +
                betas_i_delta[c] * i_covariates[j,c]
            }
            mu <- exp(log_mu);
            out <- out + dpois(data[i,j], mu, log = TRUE)*PPs[i,k]
          }
        }
      }
    } # end case of covariates on both item covariates
  }
  
  return(out)
}

# em_cycle_poisson_with_cov -------------------------------------------------------------------

em_cycle_poisson_with_cov <- function(data, item_params, weights_and_nodes,
                             p_covariates, i_covariates,
                             i_cov_on = c("alpha", "delta"),
                             fix_alphas = NULL, same_alpha = FALSE,
                             ctol_maxstep = 1e-8) {
    if (!is.null(fix_alphas)) {
      # fix alphas to the provided values
      # e step
      item_params_fixa <- c(fix_alphas, item_params)
      names(item_params_fixa) <- c(paste0("alpha", 1:ncol(data)), names(item_params))
      PPs <- estep_poisson_with_cov(
        data = data,
        item_params = item_params_fixa,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on
      )
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson_with_cov_fixalphas,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_alphas = fix_alphas,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_alpha) {
      # fit the model with estimating one same alpha for all item
      # e step
      alpha <- item_params[grepl("alpha", names(item_params))]
      item_params_samea <- c(rep(alpha, ncol(data)), item_params[-which(item_params == alpha)])
      names(item_params_samea) <- c(paste0("alpha", 1:ncol(data)), 
                                   names(item_params[-which(item_params == alpha)]))
      PPs <- estep_poisson_with_cov(
        data = data,
        item_params = item_params_samea,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on
      )
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson_with_cov_samealpha,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        control = list(xtol = ctol_maxstep)
      )$x
    } else {
      # fit a full two parameter model
      # e step
      PPs <- estep_poisson_with_cov(
        data = data,
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on
      )
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson_with_cov,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        control = list(xtol = ctol_maxstep)
      )$x
    }

  return(new_item_params)
}

# marg_ll_poisson_with_cov ---------------------------------------------------------------------

marg_ll_poisson_with_cov <- function(data, item_params, weights_and_nodes, 
                                     p_covariates, i_covariates, 
                                     i_cov_on = c("alpha", "delta"),
                                     fix_alphas = NULL, same_alphas = FALSE) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # note that deltas will be a scalar if we have item parameters on delta,
  # otherwise a vector
  if (is.null(fix_alphas)) {
    # we don't have fixed values for alpha
    if (same_alphas) {
      # we have only one alpha which is constant across items
      alpha <- item_params[grepl("alpha", names(item_params))]
      alphas <- rep(alpha, n_items) 
    } else {
      # we have an alpha for each item
      alphas <- item_params[grepl("alpha", names(item_params)) &
                              !grepl("beta", names(item_params))]
      # note that alphas will be a scalar if we have item parameters on alpha,
      # otherwise a vector
    }
  } else {
    # we have fixed values for alpha
    alphas <- fix_alphas
  }
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) { # case of person covariates
    p_covariates <- as.matrix(p_covariates)
    # function to compute integral with quadrature over
    f <- function(z, data, p_cov_data, alphas, deltas, betas_p) {
      # sum_p_cov <- sum(betas_p * as.numeric(t(p_cov_data)))
      sum_p_cov <- as.numeric(p_cov_data%*%betas_p)
      out <- 0
      for (j in 1:n_items) {
        lambda <- exp(alphas[j] * z + deltas[j] + alphas[j] * sum_p_cov)
        out <- out + (dpois(data[,j], lambda, log = TRUE))
      }
      return(exp(out))
    }
    
    marg_prob <- numeric(n_persons)
    for (i in 1:n_persons) {
      marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                             data = data[i, , drop = FALSE], 
                             p_cov_data = p_covariates[i, , drop = FALSE],
                             alphas = alphas, deltas = deltas,
                             betas_p = betas_p)
    }
    ll <- sum(log(marg_prob))
    # end case of person covariates
  } else if (!is.null(i_covariates)) { # case of item covariates
    
    if (length(i_cov_on) == 1) {
      # distinguish between the cases of where we could have the item covariates
      if (i_cov_on == "delta") { # case of covariates on delta
        # function to compute integral with quadrature over
        f <- function(z, data, i_cov_data, alphas, deltas, betas_i) {
          out <- 0
          for (j in 1:n_items) {
            sum_i_cov <- sum(betas_i * as.numeric(t(i_cov_data[j, , drop = FALSE])))
            lambda <- exp(deltas + alphas[j] * z + sum_i_cov)
            # note that deltas is just a scalar in the case of item covariates on delta
            out <- out + (dpois(data[,j], lambda, log = TRUE))
          }
          return(exp(out))
        }
        
        marg_prob <- numeric(n_persons)
        for (i in 1:n_persons) {
          marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                                 data = data[i, , drop = FALSE], 
                                 i_cov_data = i_covariates,
                                 alphas = alphas, deltas = deltas,
                                 betas_i = betas_i)
        }
        ll <- sum(log(marg_prob))
        # end case of covariates on delta
      } else if (i_cov_on == "alpha") { # case of covariates on alpha
        # function to compute integral with quadrature over
        f <- function(z, data, i_cov_data, alphas, deltas, betas_i) {
          out <- 0
          for (j in 1:n_items) {
            sum_i_cov <- z * sum(betas_i * as.numeric(t(i_cov_data[j, , drop = FALSE])))
            lambda <- exp(deltas[j] + alphas * z + sum_i_cov)
            # note that alphas is just a scalar in the case of item covariates on alpha
            out <- out + (dpois(data[,j], lambda, log = TRUE))
          }
          return(exp(out))
        }
        
        marg_prob <- numeric(n_persons)
        for (i in 1:n_persons) {
          marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                                 data = data[i, , drop = FALSE], 
                                 i_cov_data = i_covariates,
                                 alphas = alphas, deltas = deltas,
                                 betas_i = betas_i)
        }
        ll <- sum(log(marg_prob))
      }
      # end case one covariate
    } else {
      # in the poisson case, we just have the alternative that we have item 
      # covariates on both item parameters
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      
      # function to compute integral with quadrature over
      f <- function(z, data, i_cov_data, alphas, deltas, betas_i_alpha, betas_i_delta) {
        out <- 0
        for (j in 1:n_items) {
          lambda <- exp(deltas + alphas * z + 
                          z * sum(as.numeric(betas_i_alpha*i_cov_data[j,,drop=FALSE])) +
                          sum(as.numeric(betas_i_delta*i_cov_data[j,,drop=FALSE])))
          # both deltas and alphas is a scalar if we have covariates here on alpha and delta
          out <- out + (dpois(data[,j], lambda, log = TRUE))
        }
        return(exp(out))
      }
      
      marg_prob <- numeric(n_persons)
      for (i in 1:n_persons) {
        marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                               data = data[i, , drop = FALSE], 
                               i_cov_data = i_covariates,
                               alphas = alphas, deltas = deltas,
                               betas_i_alpha = betas_i_alpha,
                               betas_i_delta = betas_i_delta)
      }
      ll <- sum(log(marg_prob))
    } # end case both covariates
  } # end case of item covariates
  
  return(ll)
}

# run_em_poisson_with_cov ------------------------------------------------------------------------------


run_em_poisson_with_cov <- function(data, init_params, n_nodes, 
                           p_covariates, i_covariates,
                           i_cov_on = c("alpha", "delta"),
                           thres = Inf, prob = 0,
                           maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                           convcrit = "marglik",
                           fix_alphas = NULL, same_alpha = FALSE) {
  # i_cov_on argument: the default here is item covariates on all model parameters,
  # however, one could also choose just one of the model parameters, just given das a string
  
  # get nodes and weights for GH quadrature
  # weights_and_nodes <- gaussHermiteData(n_nodes)
  # weights_and_nodes$x <- weights_and_nodes$x * sqrt(2)
  # weights_and_nodes$w <- weights_and_nodes$w / sqrt(pi)
  weights_and_nodes<- quad_rule(n_nodes, thres = thres,prob = prob)
  
  new_params <- init_params
  conv <- FALSE
  iter <- 1
  
  new_ll <- 0
  marg_lls <- c()
  
  print("Start estimation...")
  
  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    old_params <- new_params
    new_params <- em_cycle_poisson_with_cov(
      data = data,
      item_params = old_params,
      weights_and_nodes = weights_and_nodes,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      fix_alphas = fix_alphas, 
      same_alpha = same_alpha,
      ctol_maxstep = ctol_maxstep
    )
    
    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll_poisson_with_cov(
        data = as.matrix(data),
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on, 
        fix_alphas = fix_alphas, 
        same_alphas = same_alpha)
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll_poisson_with_cov(
        data = as.matrix(data),
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        fix_alphas = fix_alphas, 
        same_alphas = same_alpha)
      marg_lls[iter] <- marg_ll
      #plot(marg_lls)
      #print(marg_lls)
    }
    
    iter <- iter + 1
  }
  
  print("Done!")
  
  # model_vcov <- compute_vcov(
  #   item_params = new_params,
  #   weights_and_nodes = weights_and_nodes, 
  #   data = data
  # )
  # 
  # se_params <- se_from_vcov(model_vcov)
  
  out <- list(
    params = new_params,
    #    se_params = se_params,
    iter = iter,
    conv = conv,
    #    vcov = model_vcov,
    marg_ll = marg_lls
  )
  return(out)
  
}

# get_start_values_poisson_with_cov -------------------------------------------------------------------

get_start_values_poisson_with_cov <- function(data, p_covariates, i_covariates, 
                                              same_alpha = FALSE, i_cov_on = c("alpha", "delta")) {
  
  # we just start with covariate weights set at 0
  if(!is.null(p_covariates)) { # for person covariates
    init_betas_p <- rep(0, ncol(p_covariates))
    init_deltas <- log(apply(data, 2, mean))
    
    if (same_alpha) {
      # just one alpha for all items
      init_alphas <- c()
      for (i in 1:ncol(data)) {
        init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
      }
      init_alphas <- mean(init_alphas)
    } else {
      # different alpha for each item
      init_alphas <- c()
      for (i in 1:ncol(data)) {
        init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
      }
    }
    
    start_values <- c(init_alphas, init_deltas, init_betas_p)
    names(start_values) <- c(
      paste0("alpha", 1:length(init_alphas)),
      paste0("delta", 1:length(init_deltas)),
      paste0("beta_p", 1:length(init_betas_p))
      )
  } else if (!is.null(i_covariates)) { # for item covariates
    # distinguish between on which item parameter we have item covariates
    if (length(i_cov_on) == 1) {
      init_betas_i <- rep(0, ncol(i_covariates))
      if (i_cov_on == "delta") {
        init_deltas <- log(mean(apply(data, 2, mean)))
        # note that for item covariates on delta, we will then have just one delta
        
        if (same_alpha) {
          # just one alpha for all items
          init_alphas <- c()
          for (i in 1:ncol(data)) {
            init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
          }
          init_alphas <- mean(init_alphas)
        } else {
          # different alpha for each item
          init_alphas <- c()
          for (i in 1:ncol(data)) {
            init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
          }
        }
      } else if (i_cov_on == "alpha") {
        init_deltas <- log(apply(data, 2, mean))
        
        # note that if have item covariates on alpha, we can't have any of the
        # constraints (because with fixed alphas their values are known anyways)
        # and with same alphas, they have to be the same and with item covariates with
        # different values for different items, they won't be the same
        
        # for covaruates on alpha, we just have one alpha
        init_alphas <- c()
        for (i in 1:ncol(data)) {
          init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
        }
        init_alphas <- mean(init_alphas)
      }
      start_values <- c(init_alphas, init_deltas, init_betas_i)
      names(start_values) <- c(
        paste0("alpha", 1:length(init_alphas)),
        paste0("delta", 1:length(init_deltas)),
        paste0("beta_i", 1:length(init_betas_i))
      )
    } else {
      # for poisson case, the alternative is just having covaraites on both
      init_betas_i <- rep(0, 2*ncol(i_covariates))
      
      # if we have item covaraites on alpha and delta, we can't have the constraint
      # on alpha that all alphas should be the same (because predicting them through
      # covariates would imply that they are different for items with different
      # covaraite values)
      
      init_deltas <- mean(log(apply(data, 2, mean)))
      
      # for covaruates on alpha, we just have one alpha
      init_alphas <- c()
      for (i in 1:ncol(data)) {
        init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
      }
      init_alphas <- mean(init_alphas)
      
      start_values <- c(init_alphas, init_deltas, init_betas_i)
      names(start_values) <- c(
        paste0("alpha", 1:length(init_alphas)),
        paste0("delta", 1:length(init_deltas)),
        paste0("beta_i_alpha", 1:(length(init_betas_i)/2)),
        paste0("beta_i_delta", 1:(length(init_betas_i)/2))
      )
    }
  }
  
  return(start_values)
}







