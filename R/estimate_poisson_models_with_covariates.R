
# estep_pois_with_cov --------------------------------------------------------------------

estep_pois_with_cov <- function(data, item_params, p_covariates, i_covariates, weights_and_nodes) {
  
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
        grad_betas_i[c] <- grad_betas_i[c] + i_covariates[j,c] * results_per_item[[j]]
      }
    }
  }
  
  out <- c(grad_alphas, grad_deltas, ifelse(is.null(i_covariates), grad_betas_p, grad_betas_i))
  return(out)
}

# grad_poisson_fixalphas --------------------------------------------------------------

grad_poisson_fixalphas <- function(item_params, PPs, weights_and_nodes, 
                                   data, fix_alphas) {
  data <- as.matrix(data)
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  grad_deltas <- numeric(length(deltas))
  
  for (j in 1:ncol(data)) {
    lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
    x_minus_lambda <- outer(data[,j], lambdas, "-")
    matrix_nodes <- matrix(
      weights_and_nodes$x,
      nrow = nrow(data), 
      ncol = length(weights_and_nodes$x),
      byrow = TRUE
    )
    x_minus_lambda_times_pp <- x_minus_lambda * PPs
    grad_deltas[j] <- sum(x_minus_lambda_times_pp)
  }
  
  return(grad_deltas)
}

# grad_poisson_fixalphas --------------------------------------------------------------

grad_poisson_samealpha <- function(item_params, PPs, weights_and_nodes, 
                                   data, fix_alphas) {
  data <- as.matrix(data)
  alphas <- rep(item_params[grepl("alpha", names(item_params))], ncol(data))
  deltas <- item_params[grepl("delta", names(item_params))]
  grad_deltas <- numeric(length(deltas))
  grad_alpha <- 0
  
  for (j in 1:ncol(data)) {
    lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
    x_minus_lambda <- outer(data[,j], lambdas, "-")
    matrix_nodes <- matrix(
      weights_and_nodes$x,
      nrow = nrow(data), 
      ncol = length(weights_and_nodes$x),
      byrow = TRUE
    )
    x_minus_lambda_times_pp <- x_minus_lambda * PPs
    grad_deltas[j] <- sum(x_minus_lambda_times_pp)
    grad_alpha <- grad_alpha + sum(matrix_nodes*x_minus_lambda_times_pp)
  }
  
  out <- c(grad_alpha, grad_deltas)
  return(out)
}

# em_cycle_poisson -------------------------------------------------------------------

em_cycle_poisson <- function(data, item_params, weights_and_nodes,
                     fix_alphas = NULL, same_alpha = FALSE,
                     ctol_maxstep = 1e-8) {
    if (!is.null(fix_alphas)) {
      # fix alphas to the provided values
      # e step
      item_params_fixa <- c(fix_alphas, item_params)
      names(item_params_fixa) <- c(paste0("alpha", 1:ncol(data)), names(item_params))
      PPs <- e_step_poisson(data, item_params_fixa, weights_and_nodes)
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson_fixalphas,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_alphas = fix_alphas,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_alpha) {
      # fit the model with estimating one same alpha for all item
      # e step
      alpha <- item_params[grepl("alpha", names(item_params))]
      item_params_samea <- c(rep(alpha, ncol(data)), item_params[-alpha])
      names(item_params_samea) <- c(paste0("alpha", 1:ncol(data)), 
                                   names(item_params[-alpha]))
      PPs <- e_step_poisson(data, item_params_samea, weights_and_nodes)
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson_samealpha,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        control = list(xtol = ctol_maxstep)
      )$x
    } else {
      # fit a full two parameter model
      # e step
      PPs <- e_step_poisson(data, item_params, weights_and_nodes)
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        control = list(xtol = ctol_maxstep)
      )$x
    }

  return(new_item_params)
}

# run_em_poisson ----------------------------------------------------------------------


run_em_poisson <- function(data, init_params, n_nodes, thres = Inf, prob = 0,
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
    new_params <- em_cycle_poisson(
      data, old_params, weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      fix_alphas = fix_alphas, same_alpha = same_alpha
    )
    
    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll2(
        as.matrix(data), new_params,
        weights_and_nodes, family = "poisson",
        fix_alphas = fix_alphas, same_alphas = same_alpha)
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll2(
        as.matrix(data), new_params,
        weights_and_nodes, family = "poisson",
        fix_alphas = fix_alphas, same_alphas = same_alpha)
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

# get_start_values_pois -----------------------------------------------------------------

get_start_values_pois <- function(data, same_alpha = FALSE) {
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
  
  start_values <- c(init_alphas, init_deltas)
  names(start_values) <- c(
    paste0("alpha", 1:length(init_alphas)),
    paste0("delta", 1:length(init_deltas))
  )
  return(start_values)
}







