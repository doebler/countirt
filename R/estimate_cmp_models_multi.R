
# e_step_multi --------------------------------------------------------------------

e_step_multi <- function(data, item_params, n_traits,
                         em_type = c("gh", "mc"),
                         weights_and_nodes = NULL, 
                         theta_samples = NULL) {
  # we don't need the alpha_constraints or disp_conctrsiants argument here; everything will automatically work
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
  disps <- exp(item_params[grepl("log_disp", names(item_params))])
  
  if (em_type == "gh") {
    
    PPs <- estep_multi_gh_cpp(
      data = as.matrix(data),
      alphas = alphas_matrix,
      deltas = deltas,
      disps = disps,
      nodes = weights_and_nodes$X,
      log_weights = weights_and_nodes$W,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_logZ_long = grid_logZ_long,
      grid_log_lambda_long = grid_log_lambda_long,
      max_mu = 200,
      min_mu = 0.001
    )
    
  } else if (em_type == "mc") {
    
    PPs <- estep_multi_mc_cpp(
      data = as.matrix(data),
      alphas = alphas_matrix,
      deltas = deltas,
      disps = disps,
      theta_samples = theta_samples,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_logZ_long = grid_logZ_long,
      grid_log_lambda_long = grid_log_lambda_long,
      max_mu = 200,
      min_mu = 0.001
    )
    
  }

  # output should be a matrix with N rows and K cols for GH 
  # and a matrix with N rows and n_sample cols for MC
  return(PPs)
}

# lasso_coord_descent --------------------------------------------------------------
# TODO equality constraints on alphas and disps implementieren, partial disp fixations implementieren
lasso_coord_descent <- function(item_params,
                                PPs,
                                data,
                                n_traits,
                                em_type = c("gh", "mc"), 
                                weights_and_nodes = NULL, 
                                theta_samples = NULL,
                                alpha_constraints = NULL,
                                disp_constraints = NULL,
                                penalize_lambda = 0,
                                max_iter = 1000,
                                ctol = 1e-3) {
  
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  if (!is.null(alpha_constraints)) {
    full_alphas <- alpha_constraints
    full_alphas[is.na(alpha_constraints)] <- alphas[is.na(alpha_constraints)]
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
  if (!is.null(disp_constraints)) {
    disps <- disp_constraints
  } else {
    disps <- exp(item_params[grepl("log_disp", names(item_params))])
  }
  
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
  
  # only check here in r whether we have GH or MC EM, we can use the same c++ functions as they just expect
  # a matrix of L columns for argument nodes and we have that for both GH and MC
  if (em_type == "gh") {
    # blockwise coordinate descent
    for (j in 1:ncol(data)) {
      conv <- FALSE
      iter <- 1
      new_params <- c(full_alphas_matrix[,j], deltas[j])
      while (!isTRUE(conv) && iter <= max_iter) {
        old_params <- new_params
        new_delta_j <- lasso_delta_update_cpp(
          alphas_j = full_alphas_matrix[,j], 
          delta_j = deltas[j],
          disp_j = disps[j],
          data_j = data[,j], 
          PPs = PPs,
          nodes = weights_and_nodes$X, 
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
        # use the previous delta_j like in sun et al.
        # TODO i don't know if it's not maybe better to use the most up to date delta in
        # the idea of cyclic coordinate descent
        full_alphas_matrix[,j] <- lasso_alpha_update_cpp(
          alphas_j = full_alphas_matrix[,j], 
          delta_j =  deltas[j], # new_delta_j,
          disp_j = disps[j],
          update_alphas_j = compute_alpha_gradient[,j],
          penalize_lambda = penalize_lambda,
          data_j = data[,j], 
          PPs = PPs,
          nodes = weights_and_nodes$X, 
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
        
        deltas[j] <- new_delta_j
        new_params <- c(full_alphas_matrix[,j], deltas[j])
        conv <- !any(abs(old_params - new_params) > ctol)
        iter <- iter + 1
      }
    }
  } else if (em_type == "mc") {
    # blockwise coordinate descent
    for (j in 1:ncol(data)) {
      conv <- FALSE
      iter <- 1
      new_params <- c(full_alphas_matrix[,j], deltas[j])
      while (!isTRUE(conv) && iter <= max_iter) {
        old_params <- new_params
        new_delta_j <- lasso_delta_update_cpp(
          alphas_j = full_alphas_matrix[,j], 
          delta_j = deltas[j],
          disp_j = disps[j],
          data_j = data[,j], 
          PPs = PPs,
          nodes = theta_samples, 
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
        # use the previous delta_j like in sun et al.
        # TODO i don't know if it's not maybe better to use the most up to date delta in
        # the idea of cyclic coordinate descent
        full_alphas_matrix[,j] <- lasso_alpha_update_cpp(
          alphas_j = full_alphas_matrix[,j], 
          delta_j = deltas[j],
          disp_j = disps[j],
          update_alphas_j = compute_alpha_gradient[,j],
          penalize_lambda = penalize_lambda,
          data_j = data[,j], 
          PPs = PPs,
          nodes = theta_samples, 
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
        
        deltas[j] <- new_delta_j
        new_params <- c(full_alphas_matrix[,j], deltas[j])
        conv <- !any(abs(old_params - new_params) > ctol)
        iter <- iter + 1
      }
    }
  }
  
  # output updated alphas and deltas
  new_alphas_deltas <- c(as.numeric(t(full_alphas_matrix)), deltas)
  names(new_alphas_deltas) <- names(item_params[!grepl("log_disp", names(item_params))])
  
  return(new_alphas_deltas)
}

# grad_nu_lasso --------------------------------------------------------------------

grad_nu_lasso <- function(log_disps,
                          other_params, 
                          PPs, 
                          data, 
                          em_type = c("gh", "mc"),
                          weights_and_nodes = NULL,
                          theta_samples = NULL) {
  # don't need to incorporate disp constriants here for now cause so far, this 
  # function only gets used when all disps will be estimated (so when no disp_constraints)
  # are in place
  
  # set up variable
  if (em_type == "gh") {
    n_traits <- ncol(weights_and_nodes$X)
  } else if (em_type == "mc") {
    n_traits <- ncol(theta_samples)
  }
  data <- as.matrix(data)
  alphas <- other_params[grepl("alpha", names(other_params))]
  alphas_matrix <- matrix(
    alphas,
    ncol = ncol(data),
    nrow = n_traits,
    byrow = TRUE
  )
  deltas <- other_params[grepl("delta", names(other_params))]
  disps <- exp(log_disps)
  
  if (em_type == "gh") {
    grads <- grad_nu_lasso_cpp(
      alphas = alphas_matrix, 
      deltas = deltas,
      disps = disps,
      data = data, 
      PPs = PPs,
      nodes = weights_and_nodes$X, 
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001)
  } else if (em_type == "mc") {
    grads <- grad_nu_lasso_cpp(
      alphas = alphas_matrix, 
      deltas = deltas,
      disps = disps,
      data = data, 
      PPs = PPs,
      nodes = theta_samples, 
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001)
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  }
  
  return(grads)
}

# grad_multi -----------------------------------------------------------------------

grad_multi <- function(item_params, 
                       PPs, 
                       data, 
                       n_traits,
                       em_type = c("gh", "mc"),
                       weights_and_nodes = NULL, theta_samples = NULL, 
                       penalize = c("none", "ridge"), 
                       # we do lasso penalty outside of
                       # this gradient because it needs its own algorithm
                       penalize_lambda = NULL,
                       alpha_constraints = NULL,
                       disp_constraints = NULL) {
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
  
  # alpha_constraints so far expect that all or a subset of alphas are fixed to certain values and the
  # reminaing freely estimated (so no equality constraints yet)
  # disp_constraints so far expect that disps are to be fixed to specific values (on original scale!)
  # (so no eqaulity constraints yet and no fixing only certain disps)
  
  # we input the alphas for the m step here through item_params so that optimization
  # method works, we then need to "supplement" the fixed values from alpha_constraints
  
  # so we also need to "supplement" dispersion values in the case of disp_constraints
  
  data <- as.matrix(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  if (!is.null(alpha_constraints)) {
    full_alphas <- alpha_constraints
    full_alphas[is.na(alpha_constraints)] <- alphas
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
  if (!is.null(disp_constraints)) {
    disps <- disp_constraints
    compute_disp_gradient <- is.na(disp_constraints)
    # estimate if NA, so not constrained to be fixed to a certain value
  } else {
    disps <- exp(item_params[grepl("log_disp", names(item_params))])
    compute_disp_gradient <- rep(TRUE, length(disps))
    # estimate all disps, we have no constraints
  }
  
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
  # thia can be passed to the C++ function and that can then check whether to compute
  # gradient or not; the C++ function can then directly output only gradients for those alphas
  # for which we want to have gradients
  
  # TODO hier weiter machen und disp_constraints implementieren, indem ich dann das
  # compute_disp_gradient argument nach c++ uebergebe
  
  # make distinction between no or ridge regularization
  if (penalize == "none") {
    # make distinction between GH and MC
    if (em_type == "gh") {
      grads <- grad_multi_gh_cpp(
        alphas = full_alphas_matrix,
        deltas = deltas,
        disps = disps,
        data = as.matrix(data),
        PPs = PPs,
        nodes = weights_and_nodes$X,
        comp_alpha_grad = compute_alpha_gradient,
        comp_disp_grad = compute_disp_gradient,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    } else if (em_type == "mc") {
      grads <- grad_multi_mc_cpp(
        alphas = full_alphas_matrix,
        deltas = deltas,
        disps = disps,
        data = as.matrix(data),
        PPs = PPs,
        theta_samples = theta_samples,
        comp_alpha_grad = compute_alpha_gradient,
        comp_disp_grad = compute_disp_gradient,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    }
  } else if (penalize == "ridge") {
    # make distinction between GH and MC
    if (em_type == "gh") {
      grads <- grad_multi_gh_ridge_cpp(
        alphas = full_alphas_matrix,
        deltas = deltas,
        disps = disps,
        data = as.matrix(data),
        PPs = PPs,
        nodes = weights_and_nodes$X,
        comp_alpha_grad = compute_alpha_gradient,
        comp_disp_grad = compute_disp_gradient,
        penalize_lambda = penalize_lambda,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    } else if (em_type == "mc") {
      grads <- grad_multi_mc_ridge_cpp(
        alphas = full_alphas_matrix,
        deltas = deltas,
        disps = disps,
        data = as.matrix(data),
        PPs = PPs,
        theta_samples = theta_samples,
        comp_alpha_grad = compute_alpha_gradient,
        comp_disp_grad = compute_disp_gradient,
        penalize_lambda = penalize_lambda,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    }
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  }
  
  return(grads)
}

# em_cycle_multi -------------------------------------------------------------------
# TODO implement equality constraints for alpha and generally constraints for disps
em_cycle_multi <- function(data, 
                           item_params, 
                           n_traits,
                           alpha_constraints = NULL, # so far only work without penalty or with ridge
                           disp_constraints = NULL,
                           em_type = c("gh", "mc"), 
                           weights_and_nodes = NULL, theta_samples = NULL,
                           penalize = c("none", "ridge", "lasso"), 
                           penalize_lambda = NULL, 
                           ctol_maxstep = 1e-8,
                           ctol_lasso = 1e-5) {
  # alpha_constraints should be a vector of the length of all the alpha parameters
  # with the same parameter names and provide the constraints for alphas
  # we start off by just assuming that we only have 0-constraints where we don't 
  # want certain items to load on certain factors
  # if there is no constraint on the parameter, expect it to have an entry of NA
  # if the value should be fixed to a certain value, then that value should be
  # given in the appropriate place in alpha_constraints
  
  # alpha_constraints so far expect that all or a subset of alphas are fixed to certain values and the
  # reminaing freely estimated (so no equality constraints yet)
  # disp_constraints so far expect that disps are to be fixed to specific values (on original scale!)
  # (so no eqaulity constraints yet and no fixing only certain disps)
  
  # neither constraints nor regularization need consideration in E step (always the same)
  PPs <- e_step_multi(
    data = data, 
    item_params = item_params, 
    n_traits = n_traits,
    em_type = em_type,
    weights_and_nodes = weights_and_nodes,
    theta_samples = theta_samples
  )
  
  # m step
  if (penalize == "lasso") {
    
    # blockwise cyclic coordinate descent for updating alpha and delta
    new_alphas_deltas <- lasso_coord_descent(
      item_params = item_params,
      PPs = PPs,
      data = data,
      n_traits = n_traits,
      em_type = em_type, 
      weights_and_nodes = weights_and_nodes, 
      theta_samples = theta_samples,
      penalize_lambda = penalize_lambda,
      disp_constraints = disp_constraints,
      alpha_constraints = alpha_constraints,
      max_iter = 1000, # TODO das hier oben als richtige argumente uebergeben
      ctol = ctol_lasso
    )
    
    if (!is.null(disp_constraints)) {
      new_log_disps <- log(disp_constraints)
    } else {
      # numerical optimization for updating disps
      log_disps <- item_params[grepl("log_disp", names(item_params))]
      new_log_disps <- nleqslv(
        x = log_disps,
        fn = grad_nu_lasso,
        other_params = new_alphas_deltas,
        PPs = PPs,
        data = data,
        em_type = em_type,
        weights_and_nodes = weights_and_nodes,
        theta_samples = theta_samples,
        control = list(xtol = ctol_maxstep)
      )$x
    }
    
    new_item_params <- c(new_alphas_deltas, new_log_disps)
    names(new_item_params) <- names(item_params)
  } else {
    # no penalization or ridge
    
    if (!is.null(alpha_constraints)) {
      # we estimate only a subset of alphas
      # we need to only input those item parameters here which need estimation
      # otherwise estimation won't work
      deltas <- item_params[grepl("delta", names(item_params))]
      alphas <- item_params[grepl("alpha", names(item_params))]
      alphas_m_step <- alphas[is.na(alpha_constraints)]
      names(alphas_m_step) <- names(alphas[is.na(alpha_constraints)])
      
      if (!is.null(disp_constraints)) {
        # don't estimate disps as they are fixed
        item_params_m_step <- c(alphas_m_step, deltas)
        names(item_params_m_step) <- c(names(alphas_m_step), names(deltas))
      } else {
        # estimate disps while not estimating all alphas
        log_disps <- item_params[grepl("log_disp", names(item_params))]
        item_params_m_step <- c(alphas_m_step, deltas, log_disps)
        names(item_params_m_step) <- c(names(alphas_m_step), names(deltas), names(log_disps))
      }
    } else if (!is.null(disp_constraints)) {
      # we estimate all alphas, but not the disps
      alphas <- item_params[grepl("alpha", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params))]
      item_params_m_step <- c(alphas, deltas)
      names(item_params_m_step) <- c(names(alphas), names(deltas))
    } else {
      # we estimate an explanatory 2pcmpm, so all items load onto all factors
      # and we estimate the dispersions
      item_params_m_step <- item_params
    }
    
    new_item_params <- nleqslv(
      x = item_params_m_step,
      fn = grad_multi,
      PPs = PPs,
      data = data,
      n_traits = n_traits,
      em_type = em_type,
      weights_and_nodes = weights_and_nodes,
      theta_samples = theta_samples,
      alpha_constraints = alpha_constraints,
      disp_constraints = disp_constraints,
      penalize = penalize,
      penalize_lambda = penalize_lambda,
      control = list(xtol = ctol_maxstep)
    )$x
    
    if (!is.null(alpha_constraints)) {
      # if we have constraints, then new_item_params as returned by new_item_params
      # is not yet the full set of item parameters but needs to be filled up
      alphas <- new_item_params[grepl("alpha", names(new_item_params))]
      deltas <- new_item_params[grepl("delta", names(new_item_params))]
      if (!is.null(disp_constraints)) {
        # dispersions are fixed to certain values
        log_disps <- log(disp_constraints)
      } else {
        # dispersions were freely estimating
        log_disps <- new_item_params[grepl("log_disp", names(new_item_params))]
      }
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
      new_item_params <- c(as.numeric(t(new_alphas_m)), deltas, log_disps)
      theta_names <- paste0("_theta", 1:n_traits)
      alpha_names <- paste0("alpha", 1:ncol(data))
      full_alpha_names <- c(outer(alpha_names, theta_names, "paste0"))
      names(new_item_params) <- c(
        full_alpha_names,
        paste0("delta", 1:ncol(data)),
        paste0("log_disp", 1:ncol(data))
      )
    } else if (!is.null(disp_constraints)) {
      # alphas have no constraints but we don't estimate dispersions, and instead have them
      # fixed to certain values
      alphas <- new_item_params[grepl("alpha", names(new_item_params))]
      deltas <- new_item_params[grepl("delta", names(new_item_params))]
      log_disps <- log(disp_constraints)
      new_item_params <- c(alphas, deltas, log_disps)
      names(new_item_params) <- names(item_params)
    }
  }
  
  return(new_item_params)
}


# marg_ll_multi -----------------------------------------------------------------------------
marg_ll_multi <- function(data, 
                          item_params, 
                          n_traits,
                          weights_and_nodes = NULL, theta_samples = NULL,
                          penalize = c("none", "ridge", "lasso"),
                          penalize_lambda = NULL,
                          em_type = c("gh", "mc")) {
  # don't need alpha_constraints or disp_constraints argument here as item_params include fixed values; 
  # this will work automatically
  
  n_items <- ncol(data)
  n_persons <- nrow(data)
  
  data <- as.matrix(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  
  alphas <- item_params[grepl("alpha", names(item_params))]
  # alphas will now have a number and a theta, so alpha1_theta1, alpha1_theta2, etc.
  # expect first all alphas (across all thetas) for item 1, and then all for item 2, etc.
  # restructure alphas into a matrix for easier combination with the
  # quadrature nodes, we have a column for each item and a row for each theta
  alphas_matrix <- matrix(
    alphas,
    ncol = ncol(data),
    nrow = n_traits,
    byrow = TRUE
  )
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  # for c++ do case disctinctions so that the input always has the expected format
  if (em_type == "gh") {
    # i can compute the unpenalized marginal likelihood and then just substract penalty here in R
    # if it is empose (see below, doesn't need to be em specific)
    ll <- marg_ll_multi_gh_cpp(data = as.matrix(data),
                               alphas = alphas_matrix,
                               deltas = deltas, 
                               disps = disps, 
                               nodes = weights_and_nodes$X,
                               log_weights = weights_and_nodes$W,
                               grid_mus = grid_mus,  
                               grid_nus = grid_nus, 
                               grid_logZ_long = grid_logZ_long,
                               grid_log_lambda_long = grid_log_lambda_long,
                               max_mu = 200,
                               min_mu = 0.001)
  } else if (em_type == "mc") {
    # i can compute the unpenalized marginal likelihood and then just substract penalty here in R
    # if it is empose (see below, doesn't need to be em specific)
    ll <- marg_ll_multi_mc_cpp(data = as.matrix(data),
                               alphas = alphas_matrix,
                               deltas = deltas, 
                               disps = disps, 
                               theta_samples = theta_samples,
                               grid_mus = grid_mus,  
                               grid_nus = grid_nus, 
                               grid_logZ_long = grid_logZ_long,
                               grid_log_lambda_long = grid_log_lambda_long,
                               max_mu = 200,
                               min_mu = 0.001)
  }
  
  # if a penalty is imposed substract it here
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


# run_em_multi ---------------------------------------------------------------------------

# TODO equality constraints on alpha and disp
# if i implement other types of constraints on disps, i need to adjust the compute_aic function
run_em_multi <- function(data, 
                         init_params, 
                         n_traits, 
                         alpha_constraints = NULL,
                         disp_constraints = NULL,
                         n_nodes = NULL, n_samples = NULL, 
                         em_type = c("gh", "mc"), 
                         fcov_prior = NULL,
                         truncate_grid = TRUE,
                         penalize = c("none", "ridge", "lasso"), 
                         penalize_lambda = NULL,
                         maxiter = 2000, convtol = 1e-5, 
                         ctol_maxstep = 1e-8, ctol_lasso = 1e-5,
                         n_samples_conv = 20, final_n_samples = 8000,
                         convcrit = "marglik") {
  
  # note that we now have n_nodes with nodes per dimension, so that total
  # number of quadrature points n_nodes^n_traits
  # we need to go down with n_nodes as we go up with n_traits
  # n_nodes is for gh em, n_samples is for mc em
  # truncate_grid is just meaningful for GH EM, it will truncate quadrature
  # grid to remove quadrature points with very low prior probability
  # fcov_prior is the prior for either type of em for latent trait cov matrix
  # you can modify expected trait correlation with that, we should just have
  # unit variance for latent traits
  
  # alpha_constraints should be a vector of the length of all the alpha parameters
  # with the same parameter names and provide the constraints for alphas
  # we start off by just assuming that we only have 0-constraints where we don't 
  # want certain items to load on certain factors
  # if there is no constraint on the parameter, expect it to have an entry of NA
  # if the value should be fixed to a certain value, then that value should be
  # given in the appropriate place in alpha_constraints
  
  # alpha_constraints so far expect that all or a subset of alphas are fixed to certain values and the
  # reminaing freely estimated (so no equality constraints yet)
  # disp_constraints so far expect that disps are to be fixed to specific values (on original scale!)
  # (so no eqaulity constraints yet and no fixing only certain disps)
  
  if (em_type == "gh") {
    # get nodes and weights for multivariate GH quadrature
    if (is.null(fcov_prior)) {
      weights_and_nodes <- init.quad(Q = n_traits, ip = n_nodes, prune = truncate_grid)
    } else {
      weights_and_nodes <- init.quad(
        Q = n_traits, 
        prior = fcov_prior,
        ip = n_nodes, 
        prune = truncate_grid)
    }
    # weights W are on log scale
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
    
    new_params <- em_cycle_multi(
      data = data,
      item_params = old_params, 
      n_traits = n_traits,
      weights_and_nodes = weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      alpha_constraints = alpha_constraints,
      disp_constraints = disp_constraints,
      em_type = em_type,
      theta_samples = theta_samples,
      penalize = penalize,
      penalize_lambda = penalize_lambda,
      ctol_lasso = ctol_lasso
    )
    #print(new_params)
    
    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll_multi(
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
      #plot(marg_lls)
      #print(marg_lls)
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
      marg_ll <- marg_ll_multi(
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
      #plot(marg_lls)
      #print(marg_lls)
    }
    
    iter <- iter + 1
  }
  
  if (em_type == "mc" & convcrit == "marglik" & conv) {
    theta_samples <- mvrnorm(final_n_samples, fcov_prior$mu, fcov_prior$sigma)
    new_params <- em_cycle_multi(
      data = data,
      item_params = new_params, 
      n_traits = n_traits,
      weights_and_nodes = weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      alpha_constraints = alpha_constraints,
      disp_constraints = disp_constraints,
      em_type = em_type,
      theta_samples = theta_samples,
      penalize = penalize,
      penalize_lambda =  penalize_lambda
    )
  }
  
  print("Done!")
  
  out <- list(
    params = new_params,
    constraints = list(
      alpha_constraints = alpha_constraints,
      disp_constraints = disp_constraints
    ),
    iter = iter,
    conv = conv,
    marg_ll = marg_lls,
    em_type = em_type
  )
  return(out)
}

# get_start_values_multi -----------------------------------------------------------------
# TODO equality constraints on alpha and disp
get_start_values_multi <- function(data, 
                                   n_traits, 
                                   alpha_constraints = NULL,
                                   disp_constraints = NULL,
                                   n_nodes = NULL, n_samples = NULL, 
                                   em_type = c("gh", "mc"), 
                                   fcov_prior = NULL, truncate_grid = TRUE,
                                   penalize = c("none", "ridge", "lasso"), 
                                   penalize_lambda = NULL,
                                   maxiter = 1000,
                                   nsim = 1000,
                                   convtol = 1e-5) {
  # alpha_constraints so far expect that all or a subset of alphas are fixed to certain values and the
  # reminaing freely estimated (so no equality constraints yet)
  # disp_constraints so far expect that disps are to be fixed to specific values (on original scale!)
  # (so no eqaulity constraints yet and no fixing only certain disps)
  
  # for CMP start values, we fit a Poisson model and get deltas and alphas from there
  init_values_pois <- get_start_values_pois_multi(
    data = data, 
    n_traits = n_traits, 
    alpha_constraints = alpha_constraints
  )
  fit_pois <- run_em_poisson_multi(
    data = data, 
    init_params = init_values_pois, 
    n_traits = n_traits, 
    n_nodes = n_nodes, 
    n_samples = n_samples, 
    em_type = em_type, 
    fcov_prior = fcov_prior,
    truncate_grid = truncate_grid,
    penalize = penalize, 
    penalize_lambda = penalize_lambda,
    maxiter = maxiter,
    alpha_constraints = alpha_constraints,
    convtol = convtol
  )
  init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
  # the outputted alphas are the whole matrix, with constrained alphas fixed at 
  # their constrained values, so e.g. fixed at 0; i thus need to pass alpha_constraints
  # to all subsequent functions and check which of these values need never be
  # estimated because they were fixed at a certain values
  init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  
  if (!is.null(disp_constraints)) {
    init_logdisps <- log(disp_constraints)
  } else {
    # start values for log nus
    init_logdisps <- c()
    sim_abilities <- mvrnorm(nsim, rep(0, n_traits), diag(rep(1, n_traits)))
    for (i in 1:ncol(data)) {
      alphas_for_item_i <- init_alphas[grepl(paste0("alpha", i, "_"), names(init_alphas))]
      mu <- exp(init_deltas[i] + sim_abilities %*% alphas_for_item_i)
      sim <- rpois(nsim, mu)
      init_logdisps[i] <- log((var(sim) / var(data[,i])))
    }
  }
  
  start_values <- c(init_alphas, init_deltas, init_logdisps)
  names(start_values) <- c(
    names(init_alphas),
    names(init_deltas),
    paste0("log_disp", 1:length(init_logdisps))
  )
  return(start_values)
}
