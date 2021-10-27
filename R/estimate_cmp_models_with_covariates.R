
# estep_cmp_with_cov ---------------------------------------------------------------------
estep_cmp_with_cov <- function(data, item_params, 
                               p_covariates, i_covariates,
                               weights_and_nodes) {
  # p_covariates is a matrix with the person covariates
  # i_covariates is a matrix with the item covariates
  
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) {
    PPs <- estep_cmp_with_icov_cpp(
      data = as.matrix(data),
      alphas = alphas,
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
  }
  
  # TODO allow for person and item covariates at the same time
  
  return(PPs)
}

# grad_cmp_with_cov ----------------------------------------------------------------------
grad_cmp_with_cov <- function(item_params, PPs, weights_and_nodes, data, 
                              p_covariates, i_covariates) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) { 
    grads <- grad_cmp_with_icov_cpp(
      alphas = alphas,
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
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  }
  
  return(grads)
}


# grad_cmp_with_cov_fixdisps----------------------------------------------------------
grad_cmp_with_cov_fixdisps <- function(item_params, PPs, weights_and_nodes, 
                                       data, fix_disps,
                                       p_covariates, i_covariates) {

  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  disps <- fix_disps
  n_items <- length(alphas)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) {
    grads <- grad_cmp_with_icov_fixdisps_cpp(
      alphas = alphas, 
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
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  # TODO allow for person and item covariates at the same time
  
  return(grads)
}

# grad_cmp_with_cov_fixalphas-------------------------------------------------------
grad_cmp_with_cov_fixalphas <- function(item_params, PPs, weights_and_nodes, 
                                        data, fix_alphas,
                                        p_covariates, i_covariates) {
  
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- length(alphas)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) {
    grads <- grad_cmp_with_icov_fixalphas_cpp(
      alphas = alphas, 
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
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  # TODO allow for person and item covariates at the same time
  
  return(grads)
}

# grad_cmp_with_cov_samedisps ---------------------------------------------------------
grad_cmp_with_cov_samedisps <- function(item_params, PPs, 
                                        weights_and_nodes, data,
                                        p_covariates, i_covariates) {
  
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(alphas)
  log_disp <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(rep(log_disp, n_items))
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) {
    grads <- grad_cmp_with_icov_samedisps_cpp(
      alphas = alphas, 
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
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  # TODO allow for person and item covariates at the same time
  
  return(grads)
}

# grad_cmp_with_cov_samealphas -----------------------------------------------------
grad_cmp_with_cov_samealphas <- function(item_params, PPs, 
                                         weights_and_nodes, data,
                                         p_covariates, i_covariates) {
  
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(deltas)
  alpha <- item_params[grepl("alpha", names(item_params))]
  alphas <- rep(alpha, n_items)
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) {
    grads <- grad_cmp_with_icov_samealphas_cpp(
      alphas = alphas, 
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
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  # TODO allow for person and item covariates at the same time
  
  return(grads)
}

# ell_cmp_with_cov -------------------------------------------------------------------
ell_cmp_with_cov <- function(item_params, PPs, weights_and_nodes, 
                             data, p_covariates, i_covariates) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) { 
    ell <- ell_cmp_with_icov_cpp(
      alphas = alphas,
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
  
  # TODO allow for person and item covariates at the same time
  
  return(ell)
}

# newem_em_cycle ---------------------------------------------------------------------
em_cycle_cmp_with_cov <- function(data, item_params, weights_and_nodes,
                                  p_covariates, i_covariates,
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
        i_covariates = i_covariates
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
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!same_disps & same_alphas) { 
      # prep for e step
      alpha <- item_params[grepl("alpha", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params))]
      n_items <- length(deltas)
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
        i_covariates = i_covariates
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
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_disps & !same_alphas) {
      # prep the parameters for the e-step
      alphas <- item_params[grepl("alpha", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params))]
      n_items <- length(deltas)
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
        i_covariates = i_covariates
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
        i_covariates = i_covariates
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
        i_covariates = i_covariates
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
                                 fix_disps = NULL, fix_alphas = NULL, 
                                 same_disps = FALSE, same_alphas = FALSE) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params))]
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
  
  if (is.null(i_covariates)) {
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
  } else if (is.null(p_covariates)) {
    ll <- marg_ll_cmp_with_icov_cpp(data = as.matrix(data),
                                   alphas = alphas,
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
  }
    
  return(ll)
}

# run_newem ----------------------------------------------------------------------
run_em_cmp_with_cov <- function(data, init_params, n_nodes, 
                                p_covariates, i_covariates, 
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
      data = data, item_params = old_params, weights_and_nodes = weights_and_nodes,
      p_covariates = p_covariates, i_covariates = i_covariates,
      ctol_maxstep = ctol_maxstep, m_method = m_method,
      fix_disps = fix_disps, fix_alphas = fix_alphas,
      same_disps = same_disps, same_alphas = same_alphas
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
        fix_disps = fix_disps, fix_alphas = fix_alphas,
        same_disps = same_disps, same_alphas = same_alphas)
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
        fix_disps = fix_disps, fix_alphas = fix_alphas,
        same_disps = same_disps, same_alphas = same_alphas)
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


# get_start_values -----------------------------------------------------------------

get_start_values <- function(data, nodes = 121, nsim = 1000,
                             init_disp_one = TRUE, same_alpha = FALSE) {
  # for CMP start values, we fit a Poisson model and get deltas and alphas from there
  if (same_alpha) {
    # just one alpha for all items
    init_values_pois <- get_start_values_pois(data, same_alpha = TRUE)
    fit_pois <- run_em_poisson(data, init_values_pois, nodes, same_alpha = TRUE)
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  } else {
    # different alpha for each item
    init_values_pois <- get_start_values_pois(data)
    fit_pois <- run_em_poisson(data, init_values_pois, nodes)
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  }
  
  init_logdisps<-c()
  sim_abilities=rnorm(nsim)
  for (i in 1:ncol(data)) {
    if (same_alpha) {
      mu <- exp(init_deltas[i] + init_alphas*sim_abilities)
    } else {
      mu <- exp(init_deltas[i] + init_alphas[i]*sim_abilities)
    }
    sim <- rpois(nsim, mu)
    init_logdisps[i] <- log((var(sim) / var(data[,i])))
  }
  
  start_values <- c(init_alphas, init_deltas, init_logdisps)
  names(start_values) <- c(
    paste0("alpha", 1:length(init_alphas)),
    paste0("delta", 1:length(init_deltas)),
    paste0("log_disp", 1:length(init_logdisps))
  )
  return(start_values)
}

