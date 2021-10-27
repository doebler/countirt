
# estep_poisson_with_cov --------------------------------------------------------------------

estep_poisson_with_cov <- function(data, item_params, p_covariates, i_covariates, weights_and_nodes) {
  
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
        as.numeric(p_covariates %*% betas_p),
        alphas[j] * weights_and_nodes$x + deltas[j],
        "+"
      ))
      PPs <- PPs + apply(lambdas, 2, function(x){dpois(data[,j], x, log = TRUE)})
    }
  } else if (is.null(p_covariates)) {
    # e step for item covariates
    i_covariates <- as.matrix(i_covariates)
    
    PPs <- matrix(
      log(weights_and_nodes$w),
      nrow = nrow(data),
      ncol = length(weights_and_nodes$x),
      byrow = TRUE
    )
    
    sum_icov <- as.numeric(i_covariates %*% betas_i)
    for (j in 1:ncol(data)) {
      lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + sum_icov[j])
      PPs <- PPs + outer(data[,j], lambdas, dpois, log = TRUE)
    }
  }
  
  PPs <- exp(PPs)
  
  PPs <- PPs / rowSums(PPs)
  
  # output should be a matrix with N rows and K cols
  return(PPs)
}

# grad_poisson -----------------------------------------------------------------------

grad_poisson_with_cov <- function(item_params, PPs, weights_and_nodes, data,
                                  p_covariates, i_covariates) {
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  grad_alphas <- numeric(length(alphas))
  grad_deltas <- numeric(length(deltas))
  
  if (is.null(i_covariates)) {
    # model with person covariates
    grad_betas_p <- numeric(length(betas_p))
    p_covariates <- as.matrix(p_covariates)
    
    results_per_item <- vector(mode = "list", length = ncol(data))
    for (j in 1:ncol(data)) {
      lambdas <- exp(outer(
        as.numeric(p_covariates %*% betas_p),
        alphas[j] * weights_and_nodes$x + deltas[j],
        "+"
      ))
      x_minus_lambda <- apply(lambdas, 2, function(x){data[,j] - x})
      matrix_nodes <- matrix(
        weights_and_nodes$x,
        nrow = nrow(data), 
        ncol = length(weights_and_nodes$x),
        byrow = TRUE
      )
      x_minus_lambda_times_pp <- x_minus_lambda * PPs
      grad_alphas[j] <- sum((matrix_nodes + sum(as.numeric(p_covariates %*% betas_p)))*
                              x_minus_lambda_times_pp)
      grad_deltas[j] <- sum(x_minus_lambda_times_pp)
      results_per_item[[j]] <- x_minus_lambda_times_pp
    }
    for (p in 1:length(betas_p)) {
      alpha_times_pcov <- outer(p_covariates[,p], alphas, "*") # output: NxM matrix
      # results_per_item is a list of length M of NxK matrices
      grad_betas_p[p] <- 0
      for (j in 1:length(ncol(data))) {
        grad_betas_p[p] <- grad_betas_p[p] + sum(alpha_times_pcov[,j] * results_per_item[[j]])
      }
    }
  } else if (is.null(p_covariates)) {
    # model with item covariates
    grad_betas_i <- numeric(length(betas_i))
    i_covariates <- as.matrix(i_covariates)
    sum_icov <- as.numeric(i_covariates %*% betas_i)
    
    results_per_item <- vector(mode = "list", length = ncol(data))
    for (j in 1:ncol(data)) {
      lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + sum_icov[j])
      x_minus_lambda <- outer(data[,j], lambdas, "-")
      matrix_nodes <- matrix(
        weights_and_nodes$x,
        nrow = nrow(data), 
        ncol = length(weights_and_nodes$x),
        byrow = TRUE
      )
      x_minus_lambda_times_pp <- x_minus_lambda * PPs
      grad_alphas[j] <- sum(matrix_nodes*x_minus_lambda_times_pp)
      grad_deltas[j] <- sum(x_minus_lambda_times_pp)
      results_per_item[[j]] <- x_minus_lambda_times_pp
    }
    for (c in 1:length(betas_i)) {
      grad_betas_i[c] <- 0
      for (j in 1:length(ncol(data))) {
        grad_betas_i[c] <- grad_betas_i[c] + sum(i_covariates[j,c] * results_per_item[[j]])
      }
    }
  }
  
  out <- c(grad_alphas, grad_deltas, ifelse(is.null(i_covariates), grad_betas_p, grad_betas_i))
  return(out)
}

# grad_poisson_with_cov_fixalphas --------------------------------------------------------------

grad_poisson_with_cov_fixalphas <- function(item_params, PPs, weights_and_nodes, 
                                   data, p_covariates, i_covariates,
                                   fix_alphas) {
  data <- as.matrix(data)
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  grad_deltas <- numeric(length(deltas))
  
  if (is.null(i_covariates)) {
    # model with person covariates
    grad_betas_p <- numeric(length(betas_p))
    p_covariates <- as.matrix(p_covariates)
    
    results_per_item <- vector(mode = "list", length = ncol(data))
    for (j in 1:ncol(data)) {
      lambdas <- exp(outer(
        as.numeric(p_covariates %*% betas_p),
        alphas[j] * weights_and_nodes$x + deltas[j],
        "+"
      ))
      x_minus_lambda <- apply(lambdas, 2, function(x){data[,j] - x})
      matrix_nodes <- matrix(
        weights_and_nodes$x,
        nrow = nrow(data), 
        ncol = length(weights_and_nodes$x),
        byrow = TRUE
      )
      x_minus_lambda_times_pp <- x_minus_lambda * PPs
      grad_deltas[j] <- sum(x_minus_lambda_times_pp)
      results_per_item[[j]] <- x_minus_lambda_times_pp
    }
    for (p in 1:length(betas_p)) {
      alpha_times_pcov <- outer(p_covariates[,p], alphas, "*") # output: NxM matrix
      # results_per_item is a list of length M of NxK matrices
      grad_betas_p[p] <- 0
      for (j in 1:length(ncol(data))) {
        grad_betas_p[p] <- grad_betas_p[p] + sum(alpha_times_pcov[,j] * results_per_item[[j]])
      }
    }
  } else if (is.null(p_covariates)) {
    # model with item covariates
    grad_betas_i <- numeric(length(betas_i))
    i_covariates <- as.matrix(i_covariates)
    sum_icov <- as.numeric(i_covariates %*% betas_i)
    
    results_per_item <- vector(mode = "list", length = ncol(data))
    for (j in 1:ncol(data)) {
      lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + sum_icov[j])
      x_minus_lambda <- outer(data[,j], lambdas, "-")
      matrix_nodes <- matrix(
        weights_and_nodes$x,
        nrow = nrow(data), 
        ncol = length(weights_and_nodes$x),
        byrow = TRUE
      )
      x_minus_lambda_times_pp <- x_minus_lambda * PPs
      grad_deltas[j] <- sum(x_minus_lambda_times_pp)
      results_per_item[[j]] <- x_minus_lambda_times_pp
    }
    for (c in 1:length(betas_i)) {
      grad_betas_i[c] <- 0
      for (j in 1:length(ncol(data))) {
        grad_betas_i[c] <- grad_betas_i[c] + sum(i_covariates[j,c] * results_per_item[[j]])
      }
    }
  }
  
  out <- c(grad_deltas, ifelse(is.null(i_covariates), grad_betas_p, grad_betas_i))
  return(out)
}

# grad_poisson_with_cov_samealpha --------------------------------------------------------------

grad_poisson_with_cov_samealpha <- function(item_params, PPs, weights_and_nodes, 
                                            data, p_covariates, i_covariates) {
  data <- as.matrix(data)
  alphas <- rep(item_params[grepl("alpha", names(item_params))], ncol(data))
  alpha <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  grad_deltas <- numeric(length(deltas))
  grad_alpha <- 0
  
  if (is.null(i_covariates)) {
    # model with person covariates
    grad_betas_p <- numeric(length(betas_p))
    p_covariates <- as.matrix(p_covariates)
    
    results_per_item <- vector(mode = "list", length = ncol(data))
    for (j in 1:ncol(data)) {
      lambdas <- exp(outer(
        as.numeric(p_covariates %*% betas_p),
        alphas[j] * weights_and_nodes$x + deltas[j],
        "+"
      ))
      x_minus_lambda <- apply(lambdas, 2, function(x){data[,j] - x})
      matrix_nodes <- matrix(
        weights_and_nodes$x,
        nrow = nrow(data), 
        ncol = length(weights_and_nodes$x),
        byrow = TRUE
      )
      x_minus_lambda_times_pp <- x_minus_lambda * PPs
      grad_alpha <- grad_alpha + sum((matrix_nodes + sum(as.numeric(p_covariates %*% betas_p)))*
                              x_minus_lambda_times_pp)
      grad_deltas[j] <- sum(x_minus_lambda_times_pp)
      results_per_item[[j]] <- x_minus_lambda_times_pp
    }
    for (p in 1:length(betas_p)) {
      alpha_times_pcov <- alpha * p_covariates[,p] # output: vector of length N
      # results_per_item is a list of length M of NxK matrices
      grad_betas_p[p] <- 0
      for (j in 1:length(ncol(data))) {
        grad_betas_p[p] <- grad_betas_p[p] + sum(alpha_times_pcov * results_per_item[[j]])
      }
    }
  } else if (is.null(p_covariates)) {
    # model with item covariates
    grad_betas_i <- numeric(length(betas_i))
    i_covariates <- as.matrix(i_covariates)
    sum_icov <- as.numeric(i_covariates %*% betas_i)
    
    results_per_item <- vector(mode = "list", length = ncol(data))
    for (j in 1:ncol(data)) {
      lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + sum_icov[j])
      x_minus_lambda <- outer(data[,j], lambdas, "-")
      matrix_nodes <- matrix(
        weights_and_nodes$x,
        nrow = nrow(data), 
        ncol = length(weights_and_nodes$x),
        byrow = TRUE
      )
      x_minus_lambda_times_pp <- x_minus_lambda * PPs
      grad_alpha <- grad_alpha + sum(matrix_nodes*x_minus_lambda_times_pp)
      grad_deltas[j] <- sum(x_minus_lambda_times_pp)
      results_per_item[[j]] <- x_minus_lambda_times_pp
    }
    for (c in 1:length(betas_i)) {
      grad_betas_i[c] <- 0
      for (j in 1:length(ncol(data))) {
        grad_betas_i[c] <- grad_betas_i[c] + sum(i_covariates[j,c] * results_per_item[[j]])
      }
    }
  }
  
  out <- c(grad_alpha, grad_deltas, ifelse(is.null(i_covariates), grad_betas_p, grad_betas_i))
  return(out)
}

# ell_cmp_with_cov ----------------------------------------------------------------------------

ell_poisson_with_cov <- function(item_params, PPs, weights_and_nodes, 
                                 data, p_covariates, i_covariates) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  N <- nrow(data)
  M <- ncol(data)
  K <- length(weights_and_nodes$x)
  I <- length(betas_i)
  P <- length(betas_p)
  
  out <- 0
  if (is.null(i_covariates)) { # we have person covariates
    for (k in 1:K) {
      for (i in 1:N) {
        for(j in 1:M) {
          log_mu <- alphas[j] * nodes[k] + deltas[j]
          for (p in 1:P) {
            log_mu <- log_mu + betas_p[p] * alphas[j] * p_covariates(i,p)
          }
          mu <- exp(log_mu);
          out <- out + dpois(data(i,j), mu, log = TRUE)*PPs(i,k)
        }
      }
    }
  } else if (is.null(p_covariates)) { # we have item covariates
    for (k in 1:K) {
      for (i in 1:N) {
        for(j in 1:M) {
          log_mu <- alphas[j] * nodes[k] + deltas[j]
          for (c in 1:I) {
            log_mu <- log_mu + betas_i[c] * i_covariates(j,c)
          }
          mu <- exp(log_mu);
          out <- out + dpois(data(i,j), mu, log = TRUE)*PPs(i,k)
        }
      }
    }
  }
  
  return(out)
}

# em_cycle_poisson_with_cov -------------------------------------------------------------------

em_cycle_poisson_with_cov <- function(data, item_params, weights_and_nodes,
                             p_covariates, i_covariates,
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
        i_covariates = i_covariates
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
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_alpha) {
      # fit the model with estimating one same alpha for all item
      # e step
      alpha <- item_params[grepl("alpha", names(item_params))]
      item_params_samea <- c(rep(alpha, ncol(data)), item_params[-alpha])
      names(item_params_samea) <- c(paste0("alpha", 1:ncol(data)), 
                                   names(item_params[-alpha]))
      PPs <- estep_poisson_with_cov(
        data = data,
        item_params = item_params_samea,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates
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
        i_covariates = i_covariates
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
        control = list(xtol = ctol_maxstep,
                       allowSingular = TRUE)
      )$x
    }

  return(new_item_params)
}

# marg_ll_poisson_with_cov ---------------------------------------------------------------------

marg_ll_poisson_with_cov <- function(data, item_params, weights_and_nodes, 
                                     p_covariates, i_covariates, 
                                     fix_alphas = NULL, same_alphas = FALSE) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  if (is.null(fix_alphas)) {
    # we don't have fixed values for alpha
    if (same_alphas) {
      # we have only one alpha which is constant across items
      alpha <- item_params[grepl("alpha", names(item_params))]
      alphas <- rep(alpha, n_items) 
    } else {
      # we have an alpha for each item
      alphas <- item_params[grepl("alpha", names(item_params))]
    }
  } else {
    # we have fixed values for alpha
    alphas <- fix_alphas
  }
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) { # case of person covariates
    # function to compute integral with quadrature over
    f <- function(z, data, p_cov_data, alphas, deltas, betas_p) {
      sum_p_cov <- betas_p * as.numeric(t(p_cov_data))
      out <- 0
      for (j in 1:n_items) {
        lambda <- exp(alphas[j] * z + deltas[j] + sum_p_cov)
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
  } else if (is.null(p_covariates)) { # case of item covariates
    # function to compute integral with quadrature over
    f <- function(z, data, i_cov_data, alphas, deltas, betas_i) {
      out <- 0
      for (j in 1:n_items) {
        sum_i_cov <- betas_i * as.numeric(t(i_cov_data[j, , drop = FALSE]))
        lambda <- exp(alphas[j] * z + deltas[j] + sum_i_cov)
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
  
  return(ll)
}

# run_em_poisson_with_cov ------------------------------------------------------------------------------


run_em_poisson_with_cov <- function(data, init_params, n_nodes, 
                           p_covariates, i_covariates,
                           thres = Inf, prob = 0,
                           maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                           convcrit = "marglik",
                           fix_alphas = NULL, same_alpha = FALSE) {
  
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

# get_start_values_poisson_with_cov -----------------------------------------------------------------

get_start_values_poisson_with_cov <- function(data, p_covariates, i_covariates, same_alpha = FALSE) {
  
  # we just start with covariate weights set at 0
  if(is.null(i_covariates)) { # for person covariates
    init_betas_p <- rep(0, ncol(p_covariates))
  } else if (is.null(p_covariates)) { # for item covariates
    init_betas_i <- rep(0, ncol(i_covariates))
  }
  
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
  
  start_values <- c(init_alphas, init_deltas, ifelse(is.null(i_covariates), init_betas_p, init_betas_i))
  names(start_values) <- c(
    paste0("alpha", 1:length(init_alphas)),
    paste0("delta", 1:length(init_deltas)),
    ifelse(is.null(p_covariates), 
           paste0("beta_i", 1:length(init_betas_i)),
           paste0("beta_p", 1:length(init_betas_p)))
  )
  return(start_values)
}







