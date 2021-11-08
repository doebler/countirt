
# estep_cmp_with_cov ---------------------------------------------------------------------
estep_cmp_with_cov <- function(data, item_params, 
                               p_covariates, i_covariates,
                               weights_and_nodes,
                               i_cov_on = c("alpha", "delta", "log_disp")) {
  # p_covariates is a matrix with the person covariates
  # i_covariates is a matrix with the item covariates
  
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  # in case of item covariates, deltas will be ust a scalar
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    PPs <- estep_cmp_with_pcov_cpp(
      data = as.matrix(data),
      alphas = alphas,
      deltas = deltas,
      disps = disps,
      betas = betas_p,
      p_cov_data = as.matrix(p_covariates),
      nodes = weights_and_nodes$x,
      weights = weights_and_nodes$w,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_logZ_long = grid_logZ_long,
      grid_log_lambda_long = grid_log_lambda_long,
      max_mu = 200,
      min_mu = 0.001
    )
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have the covaraites
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        PPs <- estep_cmp_with_icov_delta_cpp(
          data = as.matrix(data),
          alphas = alphas,
          delta = deltas,
          disps = disps,
          betas = betas_i,
          i_cov_data = as.matrix(i_covariates),
          nodes = weights_and_nodes$x,
          weights = weights_and_nodes$w,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_logZ_long = grid_logZ_long,
          grid_log_lambda_long = grid_log_lambda_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "alpha") {
        PPs <- estep_cmp_with_icov_alpha_cpp(
          data = as.matrix(data),
          alpha = alphas,
          deltas = deltas,
          disps = disps,
          betas = betas_i,
          i_cov_data = as.matrix(i_covariates),
          nodes = weights_and_nodes$x,
          weights = weights_and_nodes$w,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_logZ_long = grid_logZ_long,
          grid_log_lambda_long = grid_log_lambda_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } # TODO implement case where we have item covariates on log nu
    } # TODO implement case where we have covariates on all item parameters
  }
  
  return(PPs)
}

# grad_cmp_with_cov ----------------------------------------------------------------------
grad_cmp_with_cov <- function(item_params, PPs, weights_and_nodes, data, 
                              p_covariates, i_covariates,
                              i_cov_on = c("alpha", "delta", "log_disp")) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that for item covariates, deltas is just a scalar
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    grads <- grad_cmp_with_pcov_cpp(
      alphas = alphas,
      deltas = deltas,
      disps = disps,
      betas = betas_p,
      data = as.matrix(data),
      p_cov_data = as.matrix(p_covariates),
      PPs = PPs,
      nodes = weights_and_nodes$x,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001)
  } else if (!is.null(i_covariates)) { 
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on)) {
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_cpp(
          alphas = alphas,
          delta = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      } else if (i_cov_on == "alpha") {
        grads <- grad_cmp_with_icov_alpha_cpp(
          alpha = alphas,
          deltas = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      }
    } # TODO implement case where we have covariates on all parameters
    
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  }
  
  return(grads)
}

# grad_cmp_with_cov_fixdisps----------------------------------------------------------
grad_cmp_with_cov_fixdisps <- function(item_params, PPs, weights_and_nodes, 
                                       data, fix_disps,
                                       p_covariates, i_covariates,
                                       i_cov_on = c("alpha", "delta", "log_disp")) {

  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that deltas is a scalar if we have item covariates
  disps <- fix_disps
  n_items <- ncol(data)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    grads <- grad_cmp_with_pcov_fixdisps_cpp(
      alphas = alphas, 
      deltas = deltas, 
      disps = disps, 
      betas = betas_p,
      data = as.matrix(data),
      p_cov_data = as.matrix(p_covariates),
      PPs = PPs,
      nodes = weights_and_nodes$x,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001
    )
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameters we have the covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_fixdisps_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "alpha") {
        grads <- grad_cmp_with_icov_alpha_fixdisps_cpp(
          alpha = alphas, 
          deltas = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } 
      # note: we can't include item covariates on log nu if we fix dispersions
      # to a specific value
    } # TODO implement the case where we have item covariates on all parameters
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  return(grads)
}

# grad_cmp_with_cov_fixalphas-------------------------------------------------------
grad_cmp_with_cov_fixalphas <- function(item_params, PPs, weights_and_nodes, 
                                        data, fix_alphas,
                                        p_covariates, i_covariates,
                                        i_cov_on = c("alpha", "delta", "log_disp")) {
  
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that deltas is a scalar if we have item covariates
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- ncol(data)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    grads <- grad_cmp_with_pcov_fixalphas_cpp(
      alphas = alphas, 
      deltas = deltas, 
      disps = disps, 
      betas = betas_p,
      data = as.matrix(data),
      p_cov_data = as.matrix(p_covariates),
      PPs = PPs,
      nodes = weights_and_nodes$x,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001
    )
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have the covariates
    if (length(i_cov_on) == 1) {
      # note: if we have the constraint that alphas are fixed to specfic values
      # then we can't predict alphas through covariates
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_fixalphas_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } # TODO implement case where we have covariate on log nu
    } # TODO implement case where we have covaraites on all item parameters
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  return(grads)
}

# grad_cmp_with_cov_samedisps ---------------------------------------------------------
grad_cmp_with_cov_samedisps <- function(item_params, PPs, 
                                        weights_and_nodes, data,
                                        p_covariates, i_covariates,
                                        i_cov_on = c("alpha", "delta", "log_disp")) {
  
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that deltas is a scalar for item covariates
  n_items <- ncol(data)
  log_disp <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(rep(log_disp, n_items))
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    grads <- grad_cmp_with_pcov_samedisps_cpp(
      alphas = alphas, 
      deltas = deltas, 
      disps = disps, 
      betas = betas_p,
      data = as.matrix(data),
      p_cov_data = as.matrix(p_covariates),
      PPs = PPs,
      nodes = weights_and_nodes$x,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001
    )
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_samedisps_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          c_cov_data = as.matrix(c_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "alpha") {
        grads <- grad_cmp_with_icov_alpha_samedisps_cpp(
          alpha = alphas, 
          deltas = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          c_cov_data = as.matrix(c_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      }
      # note: we can't have covariates on log_nu if we have the constraint of
      # same disps as covaraites would have different values for the different
      # items and would then imply different disps
    } # TODO implement case where we have covariates on all item parameters
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 

  
  return(grads)
}

# grad_cmp_with_cov_samealphas -----------------------------------------------------
grad_cmp_with_cov_samealphas <- function(item_params, PPs, 
                                         weights_and_nodes, data,
                                         p_covariates, i_covariates,
                                         i_cov_on = c("alpha", "delta", "log_disp")) {
  
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that for item covariates, deltas is a scalar
  n_items <- ncol(data)
  alpha <- item_params[grepl("alpha", names(item_params))]
  alphas <- rep(alpha, n_items)
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    grads <- grad_cmp_with_pcov_samealphas_cpp(
      alphas = alphas, 
      deltas = deltas, 
      disps = disps, 
      betas = betas_p,
      data = as.matrix(data),
      p_cov_data = as.matrix(p_covariates),
      PPs = PPs,
      nodes = weights_and_nodes$x,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001
    )
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on) == 1) {
      # note: we can't have covariates on alpha if we have the constraint same_alpha
      # as we would have different values for the items on the different covariates,
      # implying different alphas
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_samealphas_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } # TODO hier noch den fall implementieren fuer i_cov_on == "log_nu"
    } # TODO den fall implementieren, dass wir kovaraiten auf allen parametern haben
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  return(grads)
}

# ell_cmp_with_cov -------------------------------------------------------------------
ell_cmp_with_cov <- function(item_params, PPs, weights_and_nodes, 
                             data, p_covariates, i_covariates,
                             i_cov_on = c("alpha", "delta", "log_disp") ) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that deltas is a scalar if we have item covariates
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    ell <- ell_cmp_with_pcov_cpp(
      alphas = alphas,
      deltas = deltas,
      disps = disps,
      betas = betas_p,
      data = as.matrix(data),
      p_cov_data = as.matrix(p_covariates),
      PPs = PPs,
      weights = weights_and_nodes$w,
      nodes = weights_and_nodes$x,
      grid_mus = grid_mus,
      grid_nus = grid_nus,
      grid_cmp_var_long = grid_cmp_var_long,
      grid_log_lambda_long = grid_log_lambda_long,
      grid_logZ_long = grid_logZ_long,
      max_mu = 200,
      min_mu = 0.001)
  } else if (!is.null(i_covariates)) { 
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        ell <- ell_cmp_with_icov_delta_cpp(
          alphas = alphas,
          delta = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          weights = weights_and_nodes$w,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      } else if (i_cov_on == "alpha") {
        ell <- ell_cmp_with_icov_alpha_cpp(
          alpha = alphas,
          deltas = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          weights = weights_and_nodes$w,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      }
    } # 
    
  }
  
  return(ell)
}

# newem_em_cycle ---------------------------------------------------------------------
em_cycle_cmp_with_cov <- function(data, item_params, weights_and_nodes,
                                  p_covariates, i_covariates,
                                  i_cov_on = c("alpha", "delta", "log_disp"),
                                  ctol_maxstep = 1e-8, m_method = "nleqslv",
                                  fix_disps = NULL, fix_alphas = NULL,
                                  same_disps = FALSE, same_alphas = FALSE) {
  
  if (is.null(fix_disps) & is.null(fix_alphas)) {
    if (!same_disps & !same_alphas) {
      # e step
      PPs <- estep_cmp_with_cov(
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
        fn = grad_cmp_with_cov,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!same_disps & same_alphas) { 
      # prep for e step
      alpha <- item_params[grepl("alpha", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params))]
      n_items <- ncol(data)
      log_disps <- item_params[grepl("log_disp", names(item_params))]
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps, betas_p, betas_i)
      names(item_params_samealph) <- c(
        paste0("alpha", 1:n_items),
        names(item_params[grepl("delta", names(item_params))]),
        names(item_params[grepl("log_disp", names(item_params))]),
        names(item_params[grepl("beta_p", names(item_params))]),
        names(item_params[grepl("beta_i", names(item_params))])
      )
      
      # e step
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_samealph,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_samealphas,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_disps & !same_alphas) {
      # prep the parameters for the e-step
      alphas <- item_params[grepl("alpha", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params))]
      n_items <- ncol(data)
      log_disp <- item_params[grepl("log_disp", names(item_params))]
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items), betas_p, betas_i)
      names(item_params_samedisp) <- c(
        names(item_params[grepl("alpha", names(item_params))]),
        names(item_params[grepl("delta", names(item_params))]),
        paste0("log_disp", 1:n_items),
        names(item_params[grepl("beta_p", names(item_params))]),
        names(item_params[grepl("beta_i", names(item_params))])
      )
      
      # e step
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_samedisp,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_samedisps,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  } else {
    if (!is.null(fix_disps)) {
      # prep for e step
      item_params_fixdisps <- c(item_params, log(fix_disps))
      names(item_params_fixdisps ) <- c(names(item_params), paste0("log_disp", 1:length(fix_disps)))
      # e step
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_fixdisps,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_fixdisps,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        fix_disps = fix_disps,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!is.null(fix_alphas)) {
      # prep for e step
      item_params_fixalphas <- c(fix_alphas, item_params)
      names(item_params_fixalphas) <- c(paste0("alpha", 1:length(fix_alphas)), names(item_params))
      # e step 
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_fixalphas,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_fixalphas,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        fix_alphas = fix_alphas,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  }
  
  return(new_item_params)
}

# marg_ll_cmp_with_cov --------------------------------------------------------------------------

marg_ll_cmp_with_cov <- function(data, item_params, weights_and_nodes, 
                                 p_covariates, i_covariates, 
                                 i_cov_on = c("alpha", "delta", "log_disp"),
                                 fix_disps = NULL, fix_alphas = NULL, 
                                 same_disps = FALSE, same_alphas = FALSE) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that deltas will be a scalar if we have item covariates
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
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    ll <- marg_ll_cmp_with_pcov_cpp(data = as.matrix(data),
                                   alphas = alphas,
                                   deltas = deltas, 
                                   disps = disps, 
                                   betas = betas_p,
                                   p_cov_data = as.matrix(p_covariates),
                                   nodes = weights_and_nodes$x,
                                   weights = weights_and_nodes$w,
                                   grid_mus = grid_mus,  
                                   grid_nus = grid_nus, 
                                   grid_logZ_long = grid_logZ_long,
                                   grid_log_lambda_long = grid_log_lambda_long,
                                   max_mu = 150,
                                   min_mu = 0.001)
  } else if (!is.null(i_covariates)) {
    # distinguish between on which parameters we have item covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        ll <- marg_ll_cmp_with_icov_delta_cpp(data = as.matrix(data),
                                        alphas = alphas,
                                        delta = deltas, 
                                        disps = disps, 
                                        betas = betas_i,
                                        i_cov_data = as.matrix(i_covariates),
                                        nodes = weights_and_nodes$x,
                                        weights = weights_and_nodes$w,
                                        grid_mus = grid_mus,  
                                        grid_nus = grid_nus, 
                                        grid_logZ_long = grid_logZ_long,
                                        grid_log_lambda_long = grid_log_lambda_long,
                                        max_mu = 150,
                                        min_mu = 0.001)
      } else if (i_cov_on == "alpha") {
        ll <- marg_ll_cmp_with_icov_alpha_cpp(data = as.matrix(data),
                                              alpha = alphas,
                                              deltas = deltas, 
                                              disps = disps, 
                                              betas = betas_i,
                                              i_cov_data = as.matrix(i_covariates),
                                              nodes = weights_and_nodes$x,
                                              weights = weights_and_nodes$w,
                                              grid_mus = grid_mus,  
                                              grid_nus = grid_nus, 
                                              grid_logZ_long = grid_logZ_long,
                                              grid_log_lambda_long = grid_log_lambda_long,
                                              max_mu = 150,
                                              min_mu = 0.001)
      } # TODO den fall einfuegen fuer kovariaten auf nu
    } # TODO den fall einfuegen dass wir auf allen item parametern kovaraiten haben
  }
    
  return(ll)
}

# run_newem ----------------------------------------------------------------------
run_em_cmp_with_cov <- function(data, init_params, n_nodes, 
                                p_covariates, i_covariates, 
                                i_cov_on = c("alpha", "delta", "log_disp"),
                                thres = Inf, prob = 0,
                                maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                                m_method = "nleqslv", convcrit = "marglik",
                                fix_disps = NULL, fix_alphas = NULL,
                                same_disps = FALSE, same_alphas = FALSE) {

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
    new_params <- em_cycle_cmp_with_cov(
      data = data, 
      item_params = old_params, 
      weights_and_nodes = weights_and_nodes,
      p_covariates = p_covariates, 
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      ctol_maxstep = ctol_maxstep,
      m_method = m_method,
      fix_disps = fix_disps, 
      fix_alphas = fix_alphas,
      same_disps = same_disps, 
      same_alphas = same_alphas
    )
    #print(new_params)

    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <-  marg_ll_cmp_with_cov(
        data = data,
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        p_covariates = p_covariates, 
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        fix_disps = fix_disps, 
        fix_alphas = fix_alphas,
        same_disps = same_disps, 
        same_alphas = same_alphas)
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll_cmp_with_cov(
        data = data,
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        p_covariates = p_covariates, 
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        fix_disps = fix_disps, 
        fix_alphas = fix_alphas,
        same_disps = same_disps, 
        same_alphas = same_alphas)
      marg_lls[iter] <- marg_ll
      plot(marg_lls)
      print(marg_lls)
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

# TODO covaraites auf log_disp implementieren

# get_start_values_cmp_with_cov ---------------------------------------------------------------------

get_start_values_cmp_with_cov <- function(data, 
                                          p_covariates,
                                          i_covariates,
                                          nodes = 121, nsim = 1000,
                                          same_alpha = FALSE,
                                          i_cov_on = c("alpha", "delta", "log_disp")) {
  # for CMP start values, we fit a Poisson model and get deltas and alphas 
  # and betas from there
  if (same_alpha) {
    # just one alpha for all items
    # note that we can't have same alpha together with item covariates on alpha
    init_values_pois <- get_start_values_poisson_with_cov(
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      same_alpha = TRUE,
      i_cov_on = "delta"
      )
    fit_pois <- run_em_poisson_with_cov(
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      init_params = init_values_pois, 
      n_nodes = nodes, 
      same_alpha = TRUE,
      i_cov_on = "delta"
      )
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  } else {
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "log_disp") {
        init_values_pois <- get_start_values_pois(
          data = data
        )
        fit_pois <- run_em_poisson(
          data = data,
          init_params = init_values_pois
        )
      } else {
        init_values_pois <- get_start_values_poisson_with_cov(
          data = data,
          p_covariates = p_covariates,
          i_covariates = i_covariates,
          i_cov_on = i_cov_on
        )
        fit_pois <- run_em_poisson_with_cov(
          data = data,
          p_covariates = p_covariates,
          i_covariates = i_covariates,
          init_params = init_values_pois, 
          n_nodes = nodes,
          i_cov_on = i_cov_on
        )
      }
    } # TODO handle case where we want item covariates on all parameters
    # different alpha for each item
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  }
  
  if (!is.null(p_covariates)) {
    # we have a model with person covariates
    init_betas_p <- fit_pois$params[grepl("beta_p", names(fit_pois$params))]
    init_logdisps<-c()
    sim_abilities=rnorm(nsim)
    for (i in 1:ncol(data)) { 
      if (same_alpha) {
        mu <- exp(init_deltas[i] + init_alphas*sim_abilities + 
                    init_alphas*sum(t(init_betas_p * t(p_covariates))))
        # here alphas is a scalar because we have the constraint same alpha here
      } else {
        mu <- exp(init_deltas[i] + init_alphas[i]*sim_abilities + 
                    init_alphas[i]*sum(t(init_betas_p * t(p_covariates))))
      }
      sim <- rpois(nsim, mu)
      init_logdisps[i] <- log((var(sim) / var(data[,i])))
    }
    
    start_values <- c(init_alphas, init_deltas, init_logdisps, init_betas_p)
    names(start_values) <- c(
      paste0("alpha", 1:length(init_alphas)),
      paste0("delta", 1:length(init_deltas)),
      paste0("log_disp", 1:length(init_logdisps)),
      paste0("beta_p", 1:length(init_betas_p))
    )
    
  } else if (!is.null(i_covariates)) {
    # we have a model with item covariates
    init_logdisps<-c()
    sim_abilities=rnorm(nsim)
    for (i in 1:ncol(data)) {
      # TODO here, I can reduce case distinctions as I pretty much do the same in all the cases
      # anyways
      # distinguish between on which parameters we have item covariates
      if (length(i_cov_on) == 1) {
        if (i_cov_on == "delta") {
          init_betas_i <- fit_pois$params[grepl("beta_i", names(fit_pois$params))]
          if (same_alpha) {
            mu <- exp(init_deltas + init_alphas*sim_abilities)
                        # + sum(t(init_betas_i * t(i_covariates))))
          } else {
            mu <- exp(init_deltas + init_alphas[i]*sim_abilities)
                        # + sum(t(init_betas_i * t(i_covariates))))
          }
        } else if (i_cov_on == "alpha") {
          init_betas_i <- fit_pois$params[grepl("beta_i", names(fit_pois$params))]
          # we can't have the constraint of same alphas if we have covaraites on
          # alpha because the covariates have different values for the different
          # items implying different alphas
          mu <- exp(init_deltas[i] + init_alphas*sim_abilities)
                     # + sim_abilities * sum(t(init_betas_i * t(i_covariates))))
        } else if (i_cov_on == "log_disp") {
          mu <- exp(init_deltas[i] + init_alphas*sim_abilities)
        } 
      } # TODO hier ein else einfuegen und den fall behandeln, dass ich auf allen
      # parametern kovariaten habe
      sim <- rpois(nsim, mu)
      init_logdisps[i] <- log((var(sim) / var(data[,i])))
    }
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "log_disp") {
        init_logdisps <- mean(init_logdisps)
        # we need one log nu and then the covariate weights on nu
        predict_log_disp_df <- data.frame(
          count_var = apply(data, 2, var)
        )
        predict_log_disp_df <- as.data.frame(cbind(predict_log_disp_df, i_covariates))
        # i_covariates is a matrix with I columnds for I covaraites and M rows for the values
        # of the M items on those I covariates
        fit_log_disp <- lm(paste0("count_var ~", 
                                  paste(colnames(i_covariates), collapse = "+" )),
                           data = predict_log_disp_df)
        init_betas_i <- log(fit_log_disp$coefficients[-1])
      }
    }
    
    start_values <- c(init_alphas, init_deltas, init_logdisps, init_betas_i)
    names(start_values) <- c(
      paste0("alpha", 1:length(init_alphas)),
      paste0("delta", 1:length(init_deltas)),
      paste0("log_disp", 1:length(init_logdisps)),
      paste0("beta_i", 1:length(init_betas_i))
    )
    
  }
  
  return(start_values)
}

