
# e_step_poisson_multi --------------------------------------------------------------------
e_step_poisson_multi <- function(data, item_params, n_traits,
                                 em_type = c("gh", "mc"),
                                 weights_and_nodes = NULL, theta_samples = NULL) {
  # we don't need the alpha_constraints argument here; everything will automatically work
  # TODO check that this still holds once i implement equality constraints for alpha
  # we don't need to consider regualrization here, posterior probs will be the same regardless
  
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  # restructure alphas into a matrix for easier combinatio with the
  # quadrature nodes, we have a column for each item and a row for each theta
  alphas_matrix <- matrix(
    alphas,
    ncol = ncol(data),
    nrow = n_traits,
    byrow = TRUE
  )
  deltas <- item_params[grepl("delta", names(item_params))]
  
  if (em_type == "gh") {
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
  } else if (em_type == "mc") {
    # for the PP in the MC-EM, we don't have the quadrature weights for the prior
    # in fact, all theta_samples have a prior weight of 1/n_samples so that they cancel out
    # (weighting is taken care of by sampling proportionally form prior to prior probability)
    PPs <- matrix(
      0, # initialize with 0 because we sum over items on log scale
      nrow = nrow(data),
      ncol = length(theta_samples[,1]),
      byrow = TRUE
    )
  }
  
  for (j in 1:ncol(data)) { # TODO see if i can rewrite with some form of apply
    if (em_type == "gh") {
      lambdas <- as.numeric(exp(weights_and_nodes$X %*% alphas_matrix[,j,drop= FALSE] +
                                  deltas[j]))
    } else if (em_type == "mc") {
      lambdas <- as.numeric(exp(theta_samples %*% alphas_matrix[,j,drop= FALSE] + deltas[j]))
    }
    PPs <- PPs + outer(data[,j], lambdas, dpois, log = TRUE)
  }
  
  PPs <- exp(PPs)
  
  PPs <- PPs / rowSums(PPs)
  
  # output should be a matrix with N rows and K cols for GH 
  # and a matrix with N rows and n_sample cols for MC
  return(PPs)
}

# grad_poisson_multi ------------------------------------------------------------------

grad_poisson_multi <- function(item_params, PPs, data, n_traits,
                               em_type = c("gh", "mc"),
                               weights_and_nodes = NULL, 
                               theta_samples = NULL, 
                               penalize = c("none", "ridge"), # we do lasso penalty outside of
                               # this gradient because it needs its own algorithm
                               penalize_lambda = NULL,
                               alpha_constraints = NULL) {
  # the PPs are the joint posterior probabilities here; we get them for
  # each dimension combination and each person, so they actually have the same
  # dimensionality as in the unidimensional case N x K, just that now, 
  # K = length(weights_and_nodes$W), so the no. of quadrature points depending
  # on number of trait dimensions and number of nodes per dimension
  
  # alpha_constraints should be a vector of the length of all the alpha parameters
  # with the same parameter names and provide the constraints for alphas
  # we start off by just assuming that we only have 0-constraints where we don't 
  # want certain items to load on certain factors
  # if there is no constraint on the parameter, expect it to have an entry of NA
  # if the value should be fixed to a certain value, then that value should be
  # given in the appropriate place in alpha_constraints
  # TODO in the future, implement also equality constraints via alpha_constraints
  # and maybe even allow delta_constraints
  
  # we input the alphas for the m step here through item_params so that optimization
  # method works, we then need to "supplement" the fixed values from alpha_constraints
  
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  if (!is.null(alpha_constraints)) {
    full_alphas <- alpha_constraints
    full_alphas[is.na(full_alphas)] <- alphas
    full_alphas_matrix <- matrix(
      full_alphas,
      ncol = ncol(data),
      nrow = n_traits,
      byrow = TRUE
    )
  } else {
    full_alphas_matrix <- matrix(
      alphas,
      ncol = ncol(data),
      nrow = n_traits,
      byrow = TRUE
    )
  }
  deltas <- item_params[grepl("delta", names(item_params))]
  
  # set-up gradient storage
  grad_deltas <- numeric(length(deltas))
  # we have an alpha for each dimension-item combination in the exploratory models
  grad_alphas_matrix <- matrix(
    NA,
    ncol = ncol(data),
    nrow = n_traits
  )
  # only write in gradient here for those alphas which aren't fixed, leave them at NA
  # otherwise and then filter them out with alpha_constraints_m later (can't just
  # filter for NAs caus gradient might break and produce NA's for estimated parameters
  # and then input and output lengths of the functions won't work anymore)
  
  # set-up a way to check for which alphas we need gradients
  if (!is.null(alpha_constraints)) {
    compute_alpha_gradient <- matrix(
      is.na(alpha_constraints),
      nrow = n_traits,
      ncol = ncol(data),
      byrow = TRUE
    )
  } else {
    compute_alpha_gradient <- matrix(
      TRUE,
      nrow = n_traits,
      ncol = ncol(data),
      byrow = TRUE
    )
  }
  
  # TODO see if i can rewrite with some form of apply
  for (j in 1:ncol(data)) {
    # for lambdas, we multiply each trait dimension with respective
    # alpha and sum across those, so that we have a vector of the length
    # of weights_and_nodes$W (number of quadrature points)
    if (em_type == "gh") { 
      lambdas <- as.numeric(exp(weights_and_nodes$X %*% full_alphas_matrix[,j,drop= FALSE] +
                                  deltas[j]))
      # lambdas is a vector of length K (= no. of quadrature points)
    } else if (em_type == "mc") {
      lambdas <- as.numeric(exp(theta_samples %*% full_alphas_matrix[,j,drop= FALSE] +
                                  deltas[j]))
      # lambdas is a vector of length n_samples
    }
    
    x_minus_lambda <- outer(data[,j], lambdas, "-") # N x K or N x n_samples matrix 
    x_minus_lambda_times_pp <- x_minus_lambda * PPs # N x K or N x n_samples matrix
    
    grad_deltas[j] <- sum(x_minus_lambda_times_pp) 
    # because we sum over the whole matrix here, we sum over all persons and 
    # all quadrature nodes for GH and all persons and all theta samples for MC
    
    if (em_type == "gh") {  
      # for alpha gradient we just want the node of the trait we are
      # taking the gradient for, so just q_k_l for specific l as in alpha_jl
      grad_alphas_matrix[,j] <- apply(
        # apply for every trait for this one vector of data, so this is like
        # a loop over traits just with apply
        weights_and_nodes$X, # matrix of quadrature points for traits
        2, # for each trait, which we put into output matrix, gradients are
        # trait and item specific
        function(x) {
          sum( # sum across all persons and nodes for that trait
            # from this matrix what we get is the multiplication with q_k_l
            # (where k is a running index, of factor l, so all quad. values of factor l)
            matrix(
              x,
              nrow = nrow(data), # length N
              ncol = length(x), # length K
              byrow = TRUE
              # for each observation (N rows) we fill the row with the vector
              # of length K (this is why K cols) of the nodes of factor l so that for each
              # person they are multiplied with the rest of the gradient
            )*
              x_minus_lambda_times_pp # N x K matrix
          )
          # for alpha constraints in confirmatory models: we are in the loop for 
          # a specific item j, so in the jth column of compute_alpha_gradient and we
          # need to check on which factors this item loads 
        }
      )
      if (penalize == "ridge") {
        grad_alphas_matrix[,j] <- grad_alphas_matrix[,j,drop=FALSE] -
          2*penalize_lambda*full_alphas_matrix[,j,drop=FALSE]
      } 
    } else if (em_type == "mc") {  
      #  this is the same code as in GH case, just swap nodes for theta samples
      # so for details on code, compare comments in GH case
      grad_alphas_matrix[,j] <- apply(
        theta_samples, # matrix with n_traits columns for each of which we need a gradient
        2,
        function(x) {
          sum(matrix(
            x,
            nrow = nrow(data), 
            ncol = length(x), 
            byrow = TRUE
            )*x_minus_lambda_times_pp
          )
        }
      )
      if (penalize == "ridge") {
        grad_alphas_matrix[,j] <- grad_alphas_matrix[,j,drop=FALSE] -
          2*penalize_lambda*full_alphas_matrix[,j,drop=FALSE]
      } 
    }
  }
  
  # prep for output, if we have constraints: here we check which alphas to output
  # gradients for and output only those so that optimizer works
  rownames(compute_alpha_gradient) <- paste0("theta", 1:nrow(compute_alpha_gradient))
  colnames(compute_alpha_gradient) <- paste0("alpha", 1:ncol(compute_alpha_gradient))
  element_names <- as.character(
    outer(colnames(compute_alpha_gradient), 
          rownames(compute_alpha_gradient),
          "paste", sep = "_"))
  compute_grad_alphas <- as.logical(t(compute_alpha_gradient))
  
  grad_alphas <- as.numeric(t(grad_alphas_matrix))[compute_grad_alphas]
  out <- c(grad_alphas, grad_deltas)
  return(out)
}

# lasso_alpha_update_poisson ---------------------------------------------------------
lasso_alpha_update_poisson <- function(alphas_j, delta_j,
                                       penalize_lambda, PPs,
                                       data_j, # a vector of responses just to item j
                                       weights_and_nodes = NULL,
                                       theta_samples = NULL, 
                                       em_type = c("gh", "mc")) {
  
  # alpha_j has all L alphas for item j, data_j is a scalar with the item intercept for item j
  # and data_j are the responses to item j for all N participants
  
  if (em_type == "gh") {
    n_traits <- ncol(weights_and_nodes$X)
  } else if (em_type == "mc") {
    n_traits <- ncol(theta_samples)
  }
   
  # TODO figure out whether i can do the whole coordinate cycling in here, i think that would probably
  # be best and the most convenient
  # then i need to make sure i update alphas_j throughout the cycle and cycle within this
  # function over the L traits and end up with new_alphas_j which i then output
  # then i can check if this gives me the results i want and if not, i can try not using most up to date
  # values right away - but i think it should be right-away updates
  # FIXME another question is probably whether i update just once or i update in coordinate descent
  # a couple of times or until convergence -> i think it's best that i handle that in the lasso_update
  # function, like how many coordinate cycles i wanna do
  
  # coordinate cycle
  new_alphas_j <- alphas_j
  
  for (l in 1:n_traits) {
    
    # first derivative (same as in grad_poisson_multi)

    if (em_type == "gh") { 
      lambdas <- as.numeric(exp(weights_and_nodes$X %*% new_alphas_j + delta_j))
      # i need to compute sum over all traits for the lambda, new_alphas_j always has the most up-to-date
      # vector of alphas for item j
      # lambdas is a vector of length K (= no. of quadrature points)
    } else if (em_type == "mc") {
      lambdas <- as.numeric(exp(theta_samples %*% new_alphas_j + delta_j))
      # lambdas is a vector of length n_samples
    }
    
    x_minus_lambda <- outer(data_j, lambdas, "-") # N x K or N x n_samples matrix 
    x_minus_lambda_times_pp <- x_minus_lambda * PPs # N x K or N x n_samples matrix
    
    # for explanation of following code, see grad_poisson_multi
    if (em_type == "gh") {  
      first_deriv_j <- sum(
        matrix(
            weights_and_nodes$X[,l],
            nrow = length(data_j), # length N
            ncol = length(weights_and_nodes$X[,l]), # length K
            byrow = TRUE)*
          x_minus_lambda_times_pp # N x K matrix
      )
    } else if (em_type == "mc") {  
      # this is the same code as in GH case, just swap nodes for theta samples
      first_deriv_j <-sum(
        matrix(
            theta_samples[,l],
            nrow = length(data_j), 
            ncol = length(theta_samples[,l]), 
            byrow = TRUE)*
          x_minus_lambda_times_pp
      )
    }
    
    # second derivative
    
    if (em_type == "gh") {  
      scnd_deriv_j <- -sum(
        matrix(
          lambdas*weights_and_nodes$X[,l]^2,
          nrow = length(data_j), # length N
          ncol = length(lambdas), # length K
          byrow = TRUE)* 
          PPs # N x K matrix
      )
    } else if (em_type == "mc") {  
      # this is the same code as in GH case, just swap nodes for theta samples
      scnd_deriv_j <- -sum(
        matrix(
          lambdas*theta_samples[,l]^2,
          nrow = length(data_j), # length N
          ncol = length(lambdas), # length K
          byrow = TRUE)* 
          PPs # N x K matrix
      )
    }
    
    # update 
    
    new_alphas_j[l] <- - soft_thresh(-scnd_deriv_j*new_alphas_j[l] + first_deriv_j, penalize_lambda) /
      scnd_deriv_j
  }
  
  return(new_alphas_j)
}

# lasso_delta_update_poisson ---------------------------------------------------------

lasso_delta_update_poisson <- function(delta_j, alphas_j,
                                       PPs,
                                       data_j,
                                       weights_and_nodes = NULL,
                                       theta_samples = NULL,
                                       em_type = c("gh", "mc")) {
  # alpha_j has all L alphas for item j, data_j is a scalar with the item intercept for item j
  # and data_j are the responses to item j for all N participants
  
  if (em_type == "gh") { 
    lambdas <- as.numeric(exp(weights_and_nodes$X %*% alphas_j + delta_j))
  } else if (em_type == "mc") {
    lambdas <- as.numeric(exp(theta_samples %*% alphas_j + delta_j))
  }
  
  x_minus_lambda <- outer(data_j, lambdas, "-") # N x K or N x n_samples matrix 
  x_minus_lambda_times_pp <- x_minus_lambda * PPs # N x K or N x n_samples matrix
  
  first_deriv_j <- sum(x_minus_lambda_times_pp) 
  # because we sum over the whole matrix here, we sum over all persons and 
  # all quadrature nodes for GH and all persons and all theta samples for MC
  
  scnd_deriv_j <- -sum(
    matrix(
      lambdas,
      nrow = length(data_j), # length N
      ncol = length(lambdas), # length K
      byrow = TRUE)* 
      PPs # N x K matrix
  )
  
  new_delta_j <- delta_j - (first_deriv_j / scnd_deriv_j)
  
  return(new_delta_j)
}

# lasso_update_poisson ---------------------------------------------------------------

lasso_coord_descent_poisson <- function(item_params,
                                 PPs,
                                 data,
                                 n_traits,
                                 em_type = c("gh", "mc"), 
                                 weights_and_nodes = NULL, 
                                 theta_samples = NULL,
                                 penalize_lambda = NULL,
                                 max_iter = 10,
                                 ctol = 1e-3) {
  # TODO implement alpha constraints for lasso penalty
  
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  alphas_matrix <- matrix(
      alphas,
      ncol = ncol(data),
      nrow = n_traits,
      byrow = TRUE
  )
  deltas <- item_params[grepl("delta", names(item_params))]
  
  # blockwise coordinate descent
  for (j in 1:ncol(data)) {
    conv <- FALSE
    iter <- 1
    new_params <- c(alphas_matrix[,j], deltas[j])
    while (!isTRUE(conv) && iter <= max_iter) {
      old_params <- new_params
      new_delta_j <- lasso_delta_update_poisson(
        delta_j = deltas[j], 
        alphas_j = alphas_matrix[,j], 
        PPs = PPs,
        data_j = data[,j], 
        weights_and_nodes = weights_and_nodes, 
        theta_samples = theta_samples, 
        em_type = em_type
        )
      # use the previous delta_j like in sun et al.
      # TODO i don't know if it's not maybe better to use the most up to date delta in
      # the idea of cyclic coordinate descent
      alphas_matrix[,j] <- lasso_alpha_update_poisson(
        alphas_j = alphas_matrix[,j], 
        delta_j = deltas[j],
        penalize_lambda = penalize_lambda,
        PPs = PPs, 
        data_j = data[,j], 
        weights_and_nodes = weights_and_nodes, 
        theta_samples = theta_samples,
        em_type = em_type
      )
      
      deltas[j] <- new_delta_j
      new_params <- c(alphas_matrix[,j], deltas[j])
      conv <- !any(abs(old_params - new_params) > ctol)
      iter <- iter + 1
    }
  }
  
  new_item_params <- c(as.numeric(t(alphas_matrix)), deltas)
  names(new_item_params) <- names(item_params)
  
  return(new_item_params)
}


# em_cycle_poisson -------------------------------------------------------------------
em_cycle_poisson_multi <- function(data, item_params, n_traits,
                                   em_type = c("gh", "mc"), 
                                   weights_and_nodes = NULL, 
                                   theta_samples = NULL,
                                   penalize = c("none", "ridge", "lasso"), 
                                   penalize_lambda = NULL, 
                                   alpha_constraints = NULL, # so far only work without penalty or with ridge
                                   ctol_maxstep = 1e-8) {
  # alpha_constraints should be a vector of the length of all the alpha parameters
  # with the same parameter names and provide the constraints for alphas
  # we start off by just assuming that we only have 0-constraints where we don't 
  # want certain items to load on certain factors
  # if there is no constraint on the parameter, expect it to have an entry of NA
  # if the value should be fixed to a certain value, then that value should be
  # given in the appropriate place in alpha_constraints
  # TODO in the future, implement also equality constraints via alpha_constraints
  # and maybe even allow delta_constraints
  
  # constraints on alpha don't need distinction here because item_params includes
  # the full set of alphas, including constraints, so alphas that aren't
  # estimated are just 0 in there and that works
  # regularization option doesn't need disctinction here (PPs are always calculated the same) 
  PPs <- e_step_poisson_multi(
    data = data, 
    item_params = item_params, 
    n_traits = n_traits,
    em_type = em_type,
    weights_and_nodes = weights_and_nodes,
    theta_samples = theta_samples
  )
  
  # m step
  if (penalize == "lasso") {
    # TODO implement alpha constraints in conjunction with lasso penalty
    
    new_item_params <- lasso_coord_descent_poisson(
      item_params = item_params,
      PPs = PPs,
      data = data,
      n_traits = n_traits,
      em_type = em_type, 
      weights_and_nodes = weights_and_nodes, 
      theta_samples = theta_samples,
      penalize_lambda = penalize_lambda,
      max_iter = 1000, # TODO das hier oben als richtige argumente uebergeben
      ctol = 1e-4
    )
    
  } else {
    if (!is.null(alpha_constraints)) {
      # we estimate a confirmatory 2pcmpm
      # we need to only input those item parameters here which need estimation
      # otherwise estimation won't work
      deltas <- item_params[grepl("delta", names(item_params))]
      alphas <- item_params[grepl("alpha", names(item_params))]
      alphas_m_step <- alphas[is.na(alpha_constraints)]
      names(alphas_m_step) <- names(alphas[is.na(alpha_constraints)])
      item_params_m_step <- c(alphas_m_step, deltas)
      names(item_params_m_step) <- c(names(alphas_m_step), names(deltas))
    } else {
      # we estimate an explanatory 2pcmpm, so all items load onto all factors
      item_params_m_step <- item_params
    }
    
    new_item_params <- nleqslv(
      x = item_params_m_step,
      fn = grad_poisson_multi,
      PPs = PPs,
      data = data,
      n_traits = n_traits,
      em_type = em_type,
      weights_and_nodes = weights_and_nodes,
      theta_samples = theta_samples,
      alpha_constraints = alpha_constraints,
      penalize = penalize,
      penalize_lambda = penalize_lambda,
      control = list(xtol = ctol_maxstep)
    )$x
    
    # TODO check that this works when i have alpha constraints 
    if (!is.null(alpha_constraints)) {
      # if we have constraints, then new_item_params as returned by new_item_params
      # is not yet the full set of item parameters but needs to be filled up
      alphas <- new_item_params[grepl("alpha", names(new_item_params))]
      deltas <- new_item_params[grepl("delta", names(new_item_params))]
      new_alphas_m <- matrix(
        alpha_constraints,
        nrow = n_traits,
        ncol = ncol(data),
        byrow = TRUE
      )
      rownames(new_alphas_m) <- paste0("theta", 1:nrow(new_alphas_m))
      colnames(new_alphas_m) <- paste0("alpha", 1:ncol(new_alphas_m))
      estimated_alphas <- names(new_item_params)[grepl("alpha", names(new_item_params))]
      estimated_alphas_l <- strsplit(estimated_alphas, "_")
      # for each of the estimated alpha, find the correct spot in the new_alphas_m matrix
      for (a in 1:length(estimated_alphas_l)) {
        new_alphas_m[estimated_alphas_l[[a]][2], estimated_alphas_l[[a]][1]] <- alphas[a]
      }
      new_item_params <- c(as.numeric(t(new_alphas_m)), deltas)
      theta_names <- paste0("_theta", 1:n_traits)
      alpha_names <- paste0("alpha", 1:ncol(data))
      full_alpha_names <- c(outer(alpha_names, theta_names, "paste0"))
      names(new_item_params) <- c(
        full_alpha_names,
        paste0("delta", 1:ncol(data))
      )
    }
  }
  
  return(new_item_params)
}

# marg_ll_poisson_multi -------------------------------------------------------------
marg_ll_poisson_multi <- function(data, item_params, n_traits,
                                  weights_and_nodes = NULL,
                                  theta_samples = NULL,
                                  penalize = c("none", "ridge", "lasso"),
                                  penalize_lambda = NULL, 
                                  em_type = c("gh", "mc")) {
  # as we still set start values for the full alpha matrix, including for those
  # alphas which we fix to certain values and don't estimate, we don't need the 
  # alpha_constraints argument here; everything will automatically work
  # TODO check that this still holds once i implement equality constraints for alpha
  # we don't need to consider here whether lambda is fixed or estimated, it is just provided
  # here via penalize_lambda and used here in computation of MLL
  
  n_items <- ncol(data)
  n_persons <- nrow(data)
  
  data <- as.matrix(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  
  alphas <- item_params[grepl("alpha", names(item_params))]
  # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  # restructure alphas into a matrix for easier combinatio with the
  # quadrature nodes, we have a column for each item and a row for each theta
  alphas_matrix <- matrix(
    alphas,
    ncol = ncol(data),
    nrow = n_traits,
    byrow = TRUE
  )
  
  marg_prob <- numeric(n_persons)
  
  if (em_type == "gh") {
    lambdas <- matrix(
      NA,
      ncol = n_items,
      nrow = length(weights_and_nodes$W)
    )
    for (j in 1:n_items) {
      lambdas[,j] <- as.numeric(exp(weights_and_nodes$X %*% alphas_matrix[,j,drop= FALSE] +
                                      deltas[j]))
    } # TODO durch apply swappen oder matrixmultiplikation
  } else if (em_type == "mc") {
    lambdas <- matrix(
      NA,
      ncol = n_items,
      nrow = length(theta_samples[,1])
    )
    for (j in 1:n_items) {
      lambdas[,j] <- as.numeric(exp(theta_samples %*% alphas_matrix[,j,drop= FALSE] + deltas[j]))
    } # TODO durch apply swappen oder matrixmultiplikation
  }
  
  for (i in 1:n_persons) {
    out <- 0
    for (j in 1:n_items) {
      out <- out + dpois(data[i,j], lambdas[,j], log = TRUE)
    }
    if (em_type == "gh") {
      # approximate integral
      marg_prob[i] <- sum(exp(out + weights_and_nodes$W)) 
    } else if (em_type == "mc") {
      # approximate integral
      marg_prob[i] <- mean(exp(out))
      # take mean here over theta samples; out should be of length n_samples
    }
  }
  
  ll <- sum(log(marg_prob)) # sum over persons
  
  if (penalize == "ridge") {
    if (!is.null(penalize_lambda)) {
      ll <- ll - ridge_penalty(alphas = alphas, lambda = penalize_lambda)
    }
  } else if (penalize == "lasso") {
    if (!is.null(penalize_lambda)) {
      ll <- ll - lasso_penalty(alphas = alphas, lambda = penalize_lambda)
    }
  }
  
  return(ll)
}

# run_em_poisson_multi ---------------------------------------------------------------
run_em_poisson_multi <- function(data, init_params, n_traits, 
                                 n_nodes = NULL, n_samples = NULL, 
                                 em_type = c("gh", "mc"), fcov_prior = NULL,
                                 truncate_grid = TRUE,
                                 penalize = c("none", "ridge", "lasso"), 
                                 penalize_lambda = NULL,
                                 maxiter = 2000, convtol = 1e-5, ctol_maxstep = 1e-8,
                                 n_samples_conv = 10, final_n_samples = 8000,
                                 convcrit = "marglik", # can also be "params"
                                 alpha_constraints = NULL) {
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
  # i cant do that for rotation as that expects orthogonal factor loadings to go in
  
  # alpha_constraints should be a vector of the length of all the alpha parameters
  # with the same parameter names and provide the constraints for alphas
  # we start off by just assuming that we only have 0-constraints where we don't 
  # want certain items to load on certain factors
  # if there is no constraint on the parameter, expect it to have an entry of NA
  # if the value should be fixed to a certain value, then that value should be
  # given in the appropriate place in alpha_constraints
  # TODO in the future, implement also equality constraints via alpha_constraints
  # and maybe even allow delta_constraints
   
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
  } else if (em_type == "mc") {
    if (is.null(fcov_prior)) {
      # multivariate standard normal prior
      fcov_prior <- list(
        mu = rep(0, n_traits),
        sigma = diag(rep(1, n_traits))
      )
    }
    # otherwise a fcov_prior has been specified and we draw samples from that for MC-EM
  }
  
  new_params <- init_params
  conv <- FALSE
  iter <- 1
  
  new_ll <- 0
  marg_lls <- c()
  
  print("Start estimation...")
  
  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    
    if (em_type == "mc") {
      # draw one set of theta samples for each em cycle
      theta_samples <- mvrnorm(n_samples, fcov_prior$mu, fcov_prior$sigma)
    } else {
      theta_samples <- NULL
    }
    
    old_params <- new_params
    
    new_params <- em_cycle_poisson_multi(
      data = data,
      item_params = old_params, 
      n_traits = n_traits,
      weights_and_nodes = weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      alpha_constraints = alpha_constraints,
      em_type = em_type,
      theta_samples = theta_samples,
      penalize = penalize,
      penalize_lambda = penalize_lambda
    )
    print(new_params)
    
    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll_poisson_multi(
        data = as.matrix(data), 
        item_params = new_params,
        n_traits = n_traits,
        weights_and_nodes = weights_and_nodes,
        em_type = em_type,
        theta_samples = theta_samples,
        penalize = penalize,
        penalize_lambda = penalize_lambda
        )
      marg_lls[iter] <- new_ll
      plot(marg_lls)
      print(marg_lls)
      if (em_type == "gh") {
        conv <- (abs(old_ll - new_ll) < convtol)
      } else if (em_type == "mc") {
        if (iter > (n_samples_conv+1)) {
          mean_marg_ll <- mean(marg_lls[(iter-n_samples_conv):iter])
          conv <- !any(abs(marg_lls[(iter-n_samples_conv):iter] - mean_marg_ll) > convtol)
        } 
      }
    } else if (convcrit == "params") {
      # TODO implement params convergence criterion for MC (i.e., check for stationarity rather than very little change)
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      # TODO move marg_ll computation to after convergence and just do it
      # once if i assess convergence on params
      marg_ll <- marg_ll_poisson_multi(
        data = as.matrix(data), 
        item_params = new_params,
        n_traits = n_traits,
        weights_and_nodes = weights_and_nodes,
        em_type = em_type,
        theta_samples = theta_samples,
        penalize = penalize,
        penalize_lambda = penalize_lambda
        )
      marg_lls[iter] <- marg_ll
      plot(marg_lls)
      print(marg_lls)
    }
    
    iter <- iter + 1
  }
  
  if (em_type == "mc" & convcrit == "marglik" & conv) {
    theta_samples <- mvrnorm(final_n_samples, fcov_prior$mu, fcov_prior$sigma)
    new_params <- em_cycle_poisson_multi(
      data = data,
      item_params = old_params, 
      n_traits = n_traits,
      weights_and_nodes = weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      alpha_constraints = alpha_constraints,
      em_type = em_type,
      theta_samples = theta_samples,
      penalize = penalize,
      penalize_lambda =  penalize_lambda
    )
  }
  
  print("Done!")
  
  out <- list(
    params = new_params,
    iter = iter,
    conv = conv,
    marg_ll = marg_lls,
    em_type = em_type
  )
  return(out)
  
}


# get_start_values_pois_multi -----------------------------------------------------------------
# TODO ueberlegen was ich mit moeglichen heywood cases in start values mache
get_start_values_pois_multi <- function(data, n_traits,  alpha_constraints = NULL) {
  # alpha_constraints should be a vector of the length of all the alpha parameters
  # with the same parameter names and provide the constraints for alphas
  # we start off by just assuming that we only have 0-constraints where we don't 
  # want certain items to load on certain factors
  # if there is no constraint on the parameter, expect it to have an entry of NA
  # if the value should be fixed to a certain value, then that value should be
  # given in the appropriate place in alpha_constraints
  # TODO in the future, implement also equality constraints via alpha_constraints
  # and maybe even allow delta_constraints
  
  init_deltas <- log(apply(data, 2, mean))
  
  M <- ncol(data)
  
  if (!is.null(alpha_constraints)) {
    # confirmatory model with constraints on alphas
     alpha_constraints_m <- t(matrix(
       alpha_constraints,
       nrow = n_traits,
       ncol = M,
       byrow = TRUE
     ))
     # we have entries everywhere in this matrix where we have a constraint
     for (i in 1:ncol(alpha_constraints_m)) {
       item_indices <- which(is.na(alpha_constraints_m[,i]))
       for (j in item_indices) {
         alpha_constraints_m[j, i] <- 
           cor(data[,j], apply(data[,item_indices[!item_indices == j]], 1, mean))
       }
     }
     init_alphas_matrix <- alpha_constraints_m
     # i am always going to also pass alpha_constraints to all functions
     # so that i can check which parameters have constraints and thus don't need
     # to be estimated but can instead just stay at their start values
  } else {
    # exploratory case where we have laodings of all items on all factors
    fa_log_data <- fa(log(test_data+0.01), nfactors = n_traits, rotate="varimax")
    loadings_fa_log_data <- as.matrix(fa_log_data$loadings)
    attributes(loadings_fa_log_data)$dimnames <- NULL
    attributes(loadings_fa_log_data)$class <- NULL
    init_alphas_matrix <- loadings_fa_log_data
    init_alphas_matrix[init_alphas_matrix < 0] <- 0
  }
  
  start_values <- c(as.numeric(init_alphas_matrix), init_deltas)
  theta_names <- paste0("_theta", 1:n_traits)
  alpha_names <- paste0("alpha", 1:ncol(data))
  full_alpha_names <- c(outer(alpha_names, theta_names, "paste0"))
  names(start_values) <- c(
    full_alpha_names,
    paste0("delta", 1:ncol(data))
  )
  return(start_values)

}







