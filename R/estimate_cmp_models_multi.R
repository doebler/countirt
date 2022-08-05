



# em_cycle_poisson -------------------------------------------------------------------

# TODO hier weitermachen und auf cmp anpassen
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


# marg_ll_multi -----------------------------------------------------------------------------
marg_ll_multi <- function(data, 
                          item_params, 
                          n_traits,
                          weights_and_nodes = NULL, theta_samples = NULL,
                          penalize = c("none", "ridge", "lasso"),
                          penalize_lambda = NULL,
                          em_type = c("gh", "mc")) {
  # don't need alpha_constraints argument here as item_params include fixed values; this will work automatically
  
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
                               max_mu = 150,
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
                               max_mu = 150,
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

# TODO disp_constraints, sowohl fixation als auch equality constraints
# TODO equality constraints on alpha
run_em_multi <- function(data, 
                         init_params, 
                         n_traits, 
                         alpha_constraints = NULL,
                         n_nodes = NULL, n_samples = NULL, 
                         em_type = c("gh", "mc"), fcov_prior = NULL,
                         truncate_grid = TRUE,
                         penalize = c("none", "ridge", "lasso"), 
                         penalize_lambda = NULL,
                         maxiter = 2000, convtol = 1e-5, ctol_maxstep = 1e-8,
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
  
  if (em_type == "gh") {
    # TODO actually use fcov_prior argument here
    # get nodes and weights for multivariate GH quadrature
    weights_and_nodes <- init.quad(Q = n_traits, ip = n_nodes, prune = truncate_grid)
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
      em_type = em_type,
      theta_samples = theta_samples,
      penalize = penalize,
      penalize_lambda = penalize_lambda
    )
    print(new_params)
    
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
      plot(marg_lls)
      print(marg_lls)
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

# get_start_values_multi -----------------------------------------------------------------
# TODO disp_constraints implementieren, also das wir die fixieren koennen und auch
# equality constraints einbauen koennen
# TODO equality constraints on alpha
get_start_values_multi <- function(data, 
                                   n_traits, 
                                   alpha_constraints = NULL,
                                   n_nodes = NULL, n_samples = NULL, 
                                   em_type = c("gh", "mc"), 
                                   fcov_prior = NULL, truncate_grid = TRUE,
                                   penalize = c("none", "ridge", "lasso"), 
                                   penalize_lambda = NULL,
                                   maxiter = 100,
                                   nsim = 1000) {
  # use maxiter here to not let the poisson model run til full convergence but just
  # enough to get some sensible start values as the poisson em will also already take
  # quite a bit of time here
  
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
    alpha_constraints = alpha_constraints
  )
  init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
  # the outputted alphas are the whole matrix, with constrained alphas fixed at 
  # their constrained values, so e.g. fixed at 0; i thus need to pass alpha_constraints
  # to all subsequent functions and check which of these values need never be
  # estimated because they were fixed at a certain values
  init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  
  # start values for log nus
  init_logdisps <- c()
  sim_abilities <- mvrnorm(nsim, rep(0, n_traits), diag(rep(1, n_traits)))
  for (i in 1:ncol(data)) {
    alphas_for_item_i <- init_alphas[grepl(paste0("alpha", i), )]
    mu <- exp(init_deltas[i] + alphas_for_item_i%*%sim_abilities)
    sim <- rpois(nsim, mu)
    init_logdisps[i] <- log((var(sim) / var(data[,i])))
  }
  
  start_values <- c(init_alphas, init_deltas, init_logdisps)
  names(start_values) <- c(
    names(init_alphas),
    names(init_deltas),
    paste0("log_disp", 1:length(init_logdisps))
  )
  return(start_values)
}
