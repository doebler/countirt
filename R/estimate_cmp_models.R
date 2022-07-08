
# newem_estep2 ---------------------------------------------------------------------
newem_estep2 <- function(data, item_params, weights_and_nodes, item_offset = NULL) {
  
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  PPs <- e_values_newem_cpp2(
    data = as.matrix(data),
    alphas = alphas,
    deltas = deltas,
    disps = disps,
    item_offset = item_offset,
    nodes = weights_and_nodes$x,
    weights = weights_and_nodes$w,
    grid_mus = grid_mus,
    grid_nus = grid_nus,
    grid_logZ_long = grid_logZ_long,
    grid_log_lambda_long = grid_log_lambda_long,
    max_mu = 200,
    min_mu = 0.001
  )
  
  return(PPs)
}

# grad_cmp_newem2 ------------------------------------------------------------------
grad_cmp_newem2 <- function(item_params, PPs, weights_and_nodes, data, item_offset = NULL) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  grads <- grad_cmp_newem_cpp2(
    alphas = alphas,
    deltas = deltas,
    disps = disps,
    item_offset = item_offset, 
    data = as.matrix(data),
    PPs = PPs,
    nodes = weights_and_nodes$x,
    grid_mus = grid_mus,
    grid_nus = grid_nus,
    grid_cmp_var_long = grid_cmp_var_long,
    grid_log_lambda_long = grid_log_lambda_long,
    grid_logZ_long = grid_logZ_long,
    max_mu = 200,
    min_mu = 0.001)
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  }
  
  #print(grads)
  return(grads)
}

# grad_cmp_fixdisps_newem ----------------------------------------------------------
grad_cmp_fixdisps_newem <- function(item_params, PPs, 
                                     weights_and_nodes, data,
                                     fix_disps, item_offset = NULL) {

  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  disps <- fix_disps
  n_items <- length(alphas)
  
  grads <- grad_cmp_fixdisps_newem_cpp(
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
    item_offset = item_offset,
    data = as.matrix(data),
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
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

# grad_cmp_fixalphas_newem -------------------------------------------------------
grad_cmp_fixalphas_newem <- function(item_params, PPs, 
                                    weights_and_nodes, data,
                                    fix_alphas, item_offset = NULL) {
  
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- length(alphas)
  
  grads <- grad_cmp_fixalphas_newem_cpp(
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
    item_offset = item_offset,
    data = as.matrix(data),
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
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

# grad_cmp_samedisps_newem ---------------------------------------------------------
grad_cmp_samedisps_newem <- function(item_params, PPs, 
                                     weights_and_nodes, data, 
                                     item_offset = NULL) {
  
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(alphas)
  log_disp <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(rep(log_disp, n_items))
  
  grads <- grad_cmp_samedisps_newem_cpp(
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
    item_offset = item_offset,
    data = as.matrix(data),
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
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

# grad_cmp_samealphas_newem -----------------------------------------------------
grad_cmp_samealphas_newem <- function(item_params, PPs, 
                                     weights_and_nodes, data, item_offset = NULL) {
  
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(deltas)
  alpha <- item_params[grepl("alpha", names(item_params))]
  alphas <- rep(alpha, n_items)
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  grads <- grad_cmp_samealphas_newem_cpp(
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
    item_offset = item_offset,
    data = as.matrix(data),
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
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

# TODO remove because i don't need it
# ell_cmp_newem -------------------------------------------------------------------
ell_cmp_newem <- function(item_params, e_values, weights_and_nodes, data) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  ell <- ell_cmp_newem_cpp(
    alphas = alphas,
    deltas = deltas,
    disps = disps,
    data = as.matrix(data),
    exp_abilities = e_values,
    grid_mus = grid_mus,
    grid_nus = grid_nus,
    grid_cmp_var_long = grid_cmp_var_long,
    grid_log_lambda_long = grid_log_lambda_long,
    grid_logZ_long = grid_logZ_long,
    max_mu = 200,
    min_mu = 0.001)
  
  return(ell)
}

# newem_em_cycle ---------------------------------------------------------------------
newem_em_cycle <- function(data, item_params, weights_and_nodes,
                           ctol_maxstep = 1e-8, m_method = "nleqslv",
                           fix_disps = NULL, fix_alphas = NULL,
                           same_disps = FALSE, same_alphas = FALSE,
                           item_offset = NULL) {
  # I handled empty item_offset arguments in run_newem, so that I can expect a vector
  # of the same length as number of items here; if i don't have any item_offsets,
  # that vector will just contain 0's (as many as we have items)

  # e_values <- newem_estep(data, item_params, weights_and_nodes)
  # 
  # new_item_params <- nleqslv(
  #   x = item_params,
  #   fn = grad_cmp_newem,
  #   e_values = e_values,
  #   weights_and_nodes = weights_and_nodes,
  #   data = data,
  #   control = list(xtol = ctol_maxstep)
  # )$x
  
  if (is.null(fix_disps) & is.null(fix_alphas)) {
    if (!same_disps & !same_alphas) {
      # e step
      PPs <- newem_estep2(data, item_params, weights_and_nodes, item_offset = item_offset)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_newem2,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!same_disps & same_alphas) {
      # prep the parameters for the e-step
      alpha <- item_params[grepl("alpha", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params))]
      n_items <- length(deltas)
      log_disps <- item_params[grepl("log_disp", names(item_params))]
      item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps)
      names(item_params_samealph) <- c(
        paste0("alpha", 1:n_items),
        names(item_params[grepl("delta", names(item_params))]),
        names(item_params[grepl("log_disp", names(item_params))])
      )
      
      # e step
      PPs <- newem_estep2(data, item_params_samealph, weights_and_nodes, item_offset = item_offset)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_samealphas_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_disps & !same_alphas) {
      # prep the parameters for the e-step
      alphas <- item_params[grepl("alpha", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params))]
      n_items <- length(deltas)
      log_disp <- item_params[grepl("log_disp", names(item_params))]
      item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items))
      names(item_params_samedisp) <- c(
        names(item_params[grepl("alpha", names(item_params))]),
        names(item_params[grepl("delta", names(item_params))]),
        paste0("log_disp", 1:n_items)
      )
      
      # e step
      PPs <- newem_estep2(data, item_params_samedisp, weights_and_nodes, item_offset = item_offset)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_samedisps_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  } else {
    if (!is.null(fix_disps)) {
      # e step
      params_estep <- c(item_params, log(fix_disps))
      names(params_estep) <- c(names(item_params), paste0("log_disp", 1:length(fix_disps)))
      PPs <- newem_estep2(data, params_estep, weights_and_nodes, item_offset = item_offset)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_fixdisps_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_disps = fix_disps,
        item_offset = item_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!is.null(fix_alphas)) {
      # e step
      params_estep <- c(fix_alphas, item_params)
      names(params_estep) <- c(paste0("alpha", 1:length(fix_alphas)), names(item_params))
      PPs <- newem_estep2(data, params_estep, weights_and_nodes, item_offset = item_offset)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_fixalphas_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_alphas = fix_alphas,
        item_offset = item_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  }
  
  return(new_item_params)
}

# run_newem ----------------------------------------------------------------------
run_newem <- function(data, init_params, n_nodes, thres = Inf, prob = 0,
                      maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                      m_method = "nleqslv", convcrit = "marglik",
                      fix_disps = NULL, fix_alphas = NULL,
                      same_disps = FALSE, same_alphas = FALSE,
                      item_offset = NULL) {

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
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, ncol(data))
  }

  print("Start estimation...")

  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    old_params <- new_params
    new_params <- newem_em_cycle(
      data = data, 
      item_params = old_params, 
      weights_and_nodes = weights_and_nodes,
      ctol_maxstep = ctol_maxstep, 
      m_method = m_method,
      fix_disps = fix_disps, 
      fix_alphas = fix_alphas,
      same_disps = same_disps, 
      same_alphas = same_alphas, 
      item_offset = item_offset
    )
    #print(new_params)

    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll2(
        data = data, 
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        family = "cmp",
        fix_disps = fix_disps, 
        fix_alphas = fix_alphas,
        same_disps = same_disps, 
        same_alphas = same_alphas,
        item_offset = item_offset
        )
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll2(
        data = data, 
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        family = "cmp",
        fix_disps = fix_disps, 
        fix_alphas = fix_alphas,
        same_disps = same_disps, 
        same_alphas = same_alphas, 
        item_offset = item_offset)
      marg_lls[iter] <- marg_ll
      #plot(marg_lls)
      #print(marg_lls)
    }

    iter <- iter + 1
  }

  print("Done!")

  out <- list(
    params = new_params,
    item_offset = item_offset,
    constraints = list(fix_alphas = fix_alphas, same_alphas = same_alphas,
                       fix_disps = fix_disps, same_disps = same_disps),
    iter = iter,
    conv = conv,
    marg_ll = marg_lls
  )
  return(out)

}


# get_start_values -----------------------------------------------------------------

get_start_values <- function(data, nodes = 121, nsim = 1000,
                             init_disp_one = TRUE, 
                             fix_disps = NULL, fix_alphas = NULL,
                             same_disps = FALSE, same_alpha = FALSE,
                             item_offset = NULL) {
  # for CMP start values, we fit a Poisson model and get deltas and alphas from there
  if (same_alpha) {
    # just one alpha for all items
    init_values_pois <- get_start_values_pois(data, same_alpha = TRUE, item_offset = item_offset)
    fit_pois <- run_em_poisson(data, init_values_pois, nodes, same_alpha = TRUE, item_offset = item_offset)
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  } else if (!is.null(fix_alphas)) {
    # we fix alphas to certain values so we don't need start values for them as
    # run_newem won't need them in the parameter vector, it will take them via
    # the fix_alphas argument
    init_values_pois <- get_start_values_pois(data, fix_alphas = fix_alphas, item_offset = item_offset)
    fit_pois <- run_em_poisson(data, init_values_pois, nodes, fix_alphas = fix_alphas, item_offset = item_offset)
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  } else {
    # different alpha for each item
    init_values_pois <- get_start_values_pois(data, item_offset = item_offset)
    fit_pois <- run_em_poisson(data, init_values_pois, nodes, item_offset = item_offset)
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  }
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, ncol(data))
  }
  
  if (!is.null(fix_disps)) {
    # we don't need start values for dispersions as they are fixed and so run_newem
    # won't need them to be part of the item params vector, it will take them via
    # the fix_disps argument
    start_values <- c(init_alphas, init_deltas)
    names(start_values) <- c(
      paste0("alpha", 1:length(init_alphas)),
      paste0("delta", 1:length(init_deltas))
    )
  } else {
    # we need start values for disps; if they are constrained to be equal, that
    # will be taken care of below
    init_logdisps<-c()
    sim_abilities=rnorm(nsim)
    for (i in 1:ncol(data)) {
      if (same_alpha) {
        mu <- exp(init_deltas[i] + init_alphas*sim_abilities + item_offset[i])
      } else {
        mu <- exp(init_deltas[i] + init_alphas[i]*sim_abilities + item_offset[i])
      }
      sim <- rpois(nsim, mu)
      init_logdisps[i] <- log((var(sim) / var(data[,i])))
    }
    # take mean for the constraint of equal dispersions
    if (same_disps) {
      init_logdisps <- mean(init_logdisps)
    }
    # for output, distinguish between cases of fixed alphas or not fixed alphas
    if (!is.null(fix_alphas)) {
      # fixed alphas, so don't output start values for them
      start_values <- c(init_deltas, init_logdisps)
      names(start_values) <- c(
        paste0("delta", 1:length(init_deltas)),
        paste0("log_disp", 1:length(init_logdisps))
      )
    } else {
      # we have at least one of each parameter type
      start_values <- c(init_alphas, init_deltas, init_logdisps)
      names(start_values) <- c(
        paste0("alpha", 1:length(init_alphas)),
        paste0("delta", 1:length(init_deltas)),
        paste0("log_disp", 1:length(init_logdisps))
      )
    }
  }

  return(start_values)
}

