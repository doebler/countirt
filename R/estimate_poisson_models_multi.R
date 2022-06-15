
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
    byrow = TRUE
  )
  deltas <- item_params[grepl("delta", names(item_params))]
  
  # the PPs are the joint posterior probabilities here; we get them for
  # each dimension combination and each person, so they actually have the same
  # dimensionality as in the unidimensional case N x K, just that now, 
  # K = length(weights_and_nodes$W), so the no. of quadrature points depending
  # on number of trait dimensions and number of nodes per dimension
  PPs <- matrix(
    weights_and_nodes$W,
    nrow = nrow(data),
    ncol = length(weights_and_nodes$W),
    byrow = TRUE
  )
  # NOTE in uni-dimensional case, we take the log of the prior weights here
  # but in the multidimensional GH function, we already ge the log weights in W
  
  for (j in 1:ncol(data)) {
    # TODO see if i can rewrite with some form of apply
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
    byrow = TRUE
  )
  deltas <- item_params[grepl("delta", names(item_params))]
  # we have an alpha for each dimension-item combination
  grad_alphas_matrix <- matrix(
    NA,
    ncol = ncol(data),
    nrow = ncol(weights_and_nodes$X)
  )
  grad_deltas <- numeric(length(deltas))
  
  # TODO see if i can rewrite with some form of apply
  for (j in 1:ncol(data)) {
    # for lambdas, we multiply each trait dimension with respective
    # alpha and sum across those, so that we have a vector of the length
    # of weights_and_nodes$W (number of quadrature points)
    lambdas <- as.numeric(exp(weights_and_nodes$X %*% alphas_matrix[,j,drop= FALSE] +
                                deltas[j]))
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
  
  grad_alphas <- as.numeric(t(grad_alphas_matrix))
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
    new_item_params <- nleqslv(
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

# marg_ll_poisson_multi -------------------------------------------------------------

marg_ll_poisson_multi <- function(data, item_params, weights_and_nodes, 
                                  fix_alphas = NULL, same_alphas = FALSE) {
  # TODO allow for constraints on alpha, think about that a little bit
  # as we have to have slightly more flexibility with this than in uni-dim
  # case so that we can do confirmatory models
  
  n_items <- ncol(data)
  n_persons <- nrow(data)
  
  data <- as.matrix(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  
  # TODO implement if-else cases for constraints on alpha
  alphas <- item_params[grepl("alpha", names(item_params))]
  # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  # restructure alphas into a matrix for easier combinatio with the
  # quadrature nodes, we have a column for each item and a row for each theta
  alphas_matrix <- matrix(
    alphas,
    ncol = ncol(data),
    nrow = ncol(weights_and_nodes$X),
    byrow = TRUE
  )
  
  # function to compute integral with quadrature over
  # f <- function(z, data, alphas_matrix, deltas) {
  #   out <- 0
  #   for (j in 1:n_items) {
  #     lambdas <- as.numeric(exp(sum(z * alphas_matrix[,j]) + deltas[j]))
  #     out <- out + (dpois(data[,j], lambdas, log = TRUE))
  #   }
  #   return(exp(out))
  # }
  #   
  # marg_prob <- numeric(n_persons)
  # for (i in 1:n_persons) {
  #   marg_prob[i] <- eval.quad(f, weights_and_nodes,
  #                            data = data[i, , drop = FALSE],
  #                            alphas_matrix = alphas_matrix, 
  #                            deltas = deltas)
  # }
  
  marg_prob <- numeric(n_persons)
  lambdas <- matrix(
    NA,
    ncol = n_items,
    nrow = length(weights_and_nodes$W)
  )
  for (j in 1:n_items) {
    lambdas[,j] <- as.numeric(exp(weights_and_nodes$X %*% alphas_matrix[,j,drop= FALSE] +
                                deltas[j]))
  } # TODO durch apply swappen oder matrixmultiplikation
  for (i in 1:n_persons) {
    out <- 0
    for (j in 1:n_items) {
      out <- out + dpois(data[i,j], lambdas[,j], log = TRUE)
    }
    marg_prob[i] <- sum(exp(out + weights_and_nodes$W))
  }
  
  ll <- sum(log(marg_prob))
  
  return(ll)
}

# run_em_poisson_multi ---------------------------------------------------------------
run_em_poisson_multi <- function(data, init_params, n_traits, 
                                 n_nodes = NULL, n_samples = NULL, 
                                 em_type = c("gh", "mc"), fcov_prior = NULL,
                                 truncate_grid = TRUE,
                                 maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                                 convcrit = "marglik",
                                 fix_alphas = NULL, same_alpha = FALSE) {
  # note that we now have n_nodes with nodes per dimension, so that total
  # number of quadrature points n_nodes^n_traits
  # we need to go down with n_nodes as we go up with n_traits
  # n_nodes is for gh em, n_samples is for mc em
  # truncate_grid is just meaningful for GH EM, it will truncate quadrature
  # grid to remove quadrature points with very low prior probability
  # fcov_prior is the prior for either type of em for latent trait cov matrix
  # you can modify expected trait correlation with that, we should just have
  # unit variance for latent traits
  # TODO when doing start values see if maybe i can get an estimate for
  # trait correlation to use for the prior here as it will help both
  # efficiency and accuracy if i do this right
   
  if (em_type == "gh") {
    # TODO actually use fcov_prior argument here
    # get nodes and weights for multivariate GH quadrature
    weights_and_nodes <- init.quad(Q = n_traits, ip = n_nodes, prune = truncate_grid)
    # NOTE the weights W are on log scale
    # default prior works for the m2ppcm
    # for output we have a list with X which is a matrix with nodes for
    # each trait (quadrature in rows with, traits in cols),
    # with n_nodes^n_traits quadrature points (so all combinations
    # across the nodes of the traits)
    # and quadrature weights W as a vector with n_nodes^n_traits entries
    # the weights are for the combination of traits (so joint)
  }
  
  new_params <- init_params
  conv <- FALSE
  iter <- 1
  
  new_ll <- 0
  marg_lls <- c()
  
  print("Start estimation...")
  
  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    old_params <- new_params
    # TODO em_type argument weiter durchreichen durch funktionen
    new_params <- em_cycle_poisson_multi(
      data, old_params, weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      fix_alphas = fix_alphas, same_alpha = same_alpha
    )
    print(new_params)
    
    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll_poisson_multi(
        as.matrix(data), new_params,
        weights_and_nodes,
        fix_alphas = fix_alphas, same_alphas = same_alpha)
      marg_lls[iter] <- new_ll
      plot(marg_lls)
      print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll_poisson_multi(
        as.matrix(data), new_params,
        weights_and_nodes,
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

get_start_values_pois <- function(data, n_traits, same_alpha = FALSE) {
  # TODO constraints wie same alpha oder fixed alpha einbauen
  init_deltas <- log(apply(data, 2, mean))
  
  if (same_alpha) {
    # TODO anpassen auf multi Fall
    # just one alpha for all items
    # init_alphas <- c()
    # for (i in 1:ncol(data)) {
    #   init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
    # }
    # init_alphas <- mean(init_alphas)
  } else {
    # different alpha for each item
    fa_log_data <- fa(log(test_data+0.01), nfactors = n_traits, rotate="varimax")
    loadings_fa_log_data <- as.matrix(fa_log_data$loadings)
    attributes(loadings_fa_log_data)$dimnames <- NULL
    attributes(loadings_fa_log_data)$class <- NULL
    init_alphas_matrix <- loadings_fa_log_data
    init_alphas_matrix[init_alphas_matrix < 0] <- 0
  }
  
  # TODO # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  # i need to transpose alpha matrix if i do as.numeric() because R works with column major order
  start_values <- c(as.numeric(init_alphas_matrix), init_deltas)
  theta_names <- paste0("_theta", 1:n_traits)
  alpha_names <- paste0("alpha", 1:ncol(data))
  c(outer(alpha_names, theta_names, paste0))
  names(start_values) <- c(
    paste0("alpha", 1:ncol(data), "_theta1"),
    paste0("alpha", 1:M, "_theta2"),
    paste0("alpha", 1:M, "_theta3"),
    paste0("delta", 1:ncol(data))
  )
  return(start_values)
}







