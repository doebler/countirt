
# e_step_poisson_multi --------------------------------------------------------------------

e_step_poisson_multi <- function(data, item_params, weights_and_nodes) {
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  # restructure alphas into a matrix for easier combinatio with the
  # quadrature nodes, we have a column for each item and a row for each theta
  alphas_matrix <- matrix(
    alphas,
    ncol = ncol(data),
    nrow = ncol(weights_and_nodes$X),
    byrow = FALSE
  )
  deltas <- item_params[grepl("delta", names(item_params))]
  
  # the PPs are the joint posterior probabilities here; we get them for
  # each dimension combination and each person, so they actually have the same
  # dimensionality as in the unidimensional case N x K, just that now, 
  # K = length(weights_and_nodes$W), so the no. of quadrature points depending
  # on number of trait dimensions and number of nodes per dimension
  PPs <- matrix(
    log(weights_and_nodes$W),
    nrow = nrow(data),
    ncol = length(weights_and_nodes$W),
    byrow = TRUE
  )
  
  for (j in 1:ncol(data)) {
    lambdas <- as.numeric(exp(weights_and_nodes$X %*% alphas_matrix[,j,drop= FALSE] +
                                deltas[j]))
    PPs <- PPs + outer(data[,j], lambdas, dpois, log = TRUE)
  }
  
  PPs <- exp(PPs)
  
  PPs <- PPs / rowSums(PPs)
  
  # output should be a matrix with N rows and K cols
  return(PPs)
}

# grad_poisson_multi ------------------------------------------------------------------

grad_poisson_multi <- function(item_params, PPs, weights_and_nodes, data) {
  # the PPs are the joint posterior probabilities here; we get them for
  # each dimension combination and each person, so they actually have the same
  # dimensionality as in the unidimensional case N x K, just that now, 
  # K = length(weights_and_nodes$W), so the no. of quadrature points depending
  # on number of trait dimensions and number of nodes per dimension
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  # restructure alphas into a matrix for easier combinatio with the
  # quadrature nodes, we have a column for each item and a row for each theta
  alphas_matrix <- matrix(
    alphas,
    ncol = ncol(data),
    nrow = ncol(weights_and_nodes$X),
    byrow = FALSE
  )
  deltas <- item_params[grepl("delta", names(item_params))]
  # we have an alpha for each dimension-item combination
  grad_alphas_matrix <- matrix(
    NA,
    ncol = ncol(data),
    nrow = ncol(weights_and_nodes$X)
  )
  grad_deltas <- numeric(length(deltas))
  
  for (j in 1:ncol(data)) {
    # for lambdas, we multiply each trait dimension with respective
    # alpha and sum across those, so that we have a vector of the length
    # of weights_and_nodes$W (number of quadrature points)
    lambdas <- exp(weights_and_nodes$X %*% alphas_matrix[,j,drop= FALSE] + deltas[j])
    # the output is a vector of length K (= no. of quadrature points)
    x_minus_lambda <- outer(data[,j], lambdas, "-")
    # matrix_nodes <- matrix(
    #   weights_and_nodes$x,
    #   nrow = nrow(data), 
    #   ncol = length(weights_and_nodes$x),
    #   byrow = TRUE
    # )
    x_minus_lambda_times_pp <- x_minus_lambda * PPs
    grad_deltas[j] <- sum(x_minus_lambda_times_pp)
    # for alpha gradient we just want the node of the trait we are
    # taking the gradient for, so just q_k_l for specific l as in alpha_jl
    grad_alphas_matrix[,j] <- apply(
      weights_and_nodes$X,
      2, 
      function(x) {sum(
        matrix(
          x,
          nrow = nrow(data), 
          ncol = length(x),
          byrow = TRUE
        )*
        x_minus_lambda_times_pp
        )}
    )
  }
  
  grad_alphas <- grad_alphas_matrix
  out <- c(grad_alphas, grad_deltas)
  return(out)
}

# grad_poisson_fixalphas --------------------------------------------------------------
# TODO
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

# grad_poisson_samealpha --------------------------------------------------------------
# TODO
grad_poisson_samealpha <- function(item_params, PPs, weights_and_nodes, 
                                   data) {
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

em_cycle_poisson_multi <- function(
    data, item_params, weights_and_nodes,
    fix_alphas = NULL, same_alpha = FALSE,
    ctol_maxstep = 1e-8) {
  
  # TODO
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
    # TODO
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
    # fit a full multidimensional two parameter model
    # e step
    PPs <- e_step_poisson_multi(data, item_params, weights_and_nodes)
    
    # m step
    new_item_params_multi <- nleqslv(
      x = item_params,
      fn = grad_poisson_multi,
      PPs = PPs,
      weights_and_nodes = weights_and_nodes,
      data = data,
      control = list(xtol = ctol_maxstep)
    )$x
  }
  
  return(new_item_params)
}

# run_em_poisson_multi ---------------------------------------------------------------
run_em_poisson_multi <- function(data, init_params, n_traits, n_nodes, 
                                 thres = Inf, prob = 0,
                                 maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                                 convcrit = "marglik",
                                 fix_alphas = NULL, same_alpha = FALSE) {
  # note that we now have n_nodes with nodes per dimension, so that total
  # number of quadrature points n_nodes^n_traits
  # we need to go down with n_nodes as we go up with n_traits
  
  # get nodes and weights for multivariate GH quadrature
  weights_and_nodes <- init.quad(Q = n_traits, ip = n_nodes)
  # default prior works for the m2ppcm
  # for output we have a list with X which is a matrix with nodes for
  # each trait (quadrature in rows with, traits in cols),
  # with n_nodes^n_traits quadrature points (so all combinations
  # across the nodes of the traits)
  # and quadrature weights W as a vector with n_nodes^n_traits entries
  # the weights are for the combination of traits (so joint)
  
  new_params <- init_params
  conv <- FALSE
  iter <- 1
  
  new_ll <- 0
  marg_lls <- c()
  
  print("Start estimation...")
  
  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    old_params <- new_params
    new_params <- em_cycle_poisson_multi(
      data, old_params, weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      fix_alphas = fix_alphas, same_alpha = same_alpha
    )
    
    # check for convergence
    # TODO die marg_ll2 funktion anpassen fuer multi
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
  
  out <- list(
    params = new_params,
    iter = iter,
    conv = conv,
    marg_ll = marg_lls
  )
  return(out)
  
}

# TODO startwerte fuer multi
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
  
  # TODO # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  start_values <- c(init_alphas, init_deltas)
  names(start_values) <- c(
    paste0("alpha", 1:length(init_alphas)),
    paste0("delta", 1:length(init_deltas))
  )
  return(start_values)
}







