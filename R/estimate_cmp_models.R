# newem_estep <- function(data, item_params, weights_and_nodes) {
#   # prep item parameters
#   alphas <- item_params[grepl("alpha", names(item_params))]
#   deltas <- item_params[grepl("delta", names(item_params))]
#   log_disps <- item_params[grepl("log_disp", names(item_params))]
#   disps <- exp(log_disps)
# 
#   exp_abilities <- e_values_newem_cpp(
#     data = as.matrix(data),
#     alphas = alphas,
#     deltas = deltas,
#     disps = disps,
#     nodes = weights_and_nodes$x,
#     weights = weights_and_nodes$w,
#     grid_mus = grid_mus,
#     grid_nus = grid_nus,
#     grid_logZ_long = grid_logZ_long,
#     grid_log_lambda_long = grid_log_lambda_long,
#     max_mu = 200,
#     min_mu = 0.001
#   )
# 
#   print(exp_abilities)
#   return(exp_abilities)
# }

newem_estep2 <- function(data, item_params, weights_and_nodes) {
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


# grad_cmp_newem <- function(item_params, e_values, weights_and_nodes, data) {
#   # prep item parameters
#   alphas <- item_params[grepl("alpha", names(item_params))]
#   deltas <- item_params[grepl("delta", names(item_params))]
#   log_disps <- item_params[grepl("log_disp", names(item_params))]
#   disps <- exp(log_disps)
# 
#   grads <- grad_cmp_newem_cpp(
#     alphas = alphas,
#     deltas = deltas,
#     disps = disps,
#     data = as.matrix(data),
#     exp_abilities = e_values,
#     grid_mus = grid_mus,
#     grid_nus = grid_nus,
#     grid_cmp_var_long = grid_cmp_var_long,
#     grid_log_lambda_long = grid_log_lambda_long,
#     grid_logZ_long = grid_logZ_long,
#     max_mu = 200,
#     min_mu = 0.001)
# 
#   if (any(is.na(grads))) {
#     stop("Gradient contained NA", paste0(grads, collapse = ","),
#          paste0(item_params, collapse = ","))
#   }
# 
#   print(grads)
#   return(grads)
# }

grad_cmp_newem2 <- function(item_params, PPs, weights_and_nodes, data) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  grads <- grad_cmp_newem_cpp2(
    alphas = alphas,
    deltas = deltas,
    disps = disps,
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

grad_cmp_fixdisps_newem <- function(item_params, PPs, 
                                     weights_and_nodes, data,
                                     fix_disps) {

  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  disps <- fix_disps
  n_items <- length(alphas)
  
  grads <- grad_cmp_fixdisps_newem_cpp(
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
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

grad_cmp_fixalphas_newem <- function(item_params, PPs, 
                                    weights_and_nodes, data,
                                    fix_alphas) {
  
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- length(alphas)
  
  grads <- grad_cmp_fixalphas_newem_cpp(
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
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

grad_cmp_samedisps_newem <- function(item_params, PPs, 
                                     weights_and_nodes, data) {
  
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(alphas)
  log_disp <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(rep(log_disp, n_items))
  
  grads <- grad_cmp_samedisps_newem_cpp(
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
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

grad_cmp_samealphas_newem <- function(item_params, PPs, 
                                     weights_and_nodes, data) {
  
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

newem_em_cycle <- function(data, item_params, weights_and_nodes,
                           ctol_maxstep = 1e-8, m_method = "nleqslv",
                           fix_disps = NULL, fix_alphas = NULL,
                           same_disps = FALSE, same_alphas = FALSE) {

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
      PPs <- newem_estep2(data, item_params, weights_and_nodes)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_newem2,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
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
      PPs <- newem_estep2(data, item_params_samealph, weights_and_nodes)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_samealphas_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
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
      PPs <- newem_estep2(data, item_params_samedisp, weights_and_nodes)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_samedisps_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  } else {
    if (!is.null(fix_disps)) {
      # e step
      params_estep <- c(item_params, log(fix_disps))
      names(params_estep) <- c(names(item_params), paste0("log_disp", 1:length(fix_disps)))
      PPs <- newem_estep2(data, params_estep, weights_and_nodes)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_fixdisps_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_disps = fix_disps,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!is.null(fix_alphas)) {
      # e step
      params_estep <- c(fix_alphas, item_params)
      names(params_estep) <- c(paste0("alpha", 1:length(fix_alphas)), names(item_params))
      PPs <- newem_estep2(data, params_estep, weights_and_nodes)
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_fixalphas_newem,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_alphas = fix_alphas,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  }
  
  return(new_item_params)
}

run_newem <- function(data, init_params, n_nodes, thres = Inf, prob = 0,
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
    new_params <- newem_em_cycle(
      data, old_params, weights_and_nodes,
      ctol_maxstep = ctol_maxstep, m_method = m_method,
      fix_disps = fix_disps, fix_alphas = fix_alphas,
      same_disps = same_disps, same_alphas = same_alphas
    )
    #print(new_params)

    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll2(
        data, new_params,
        weights_and_nodes, family = "cmp",
        fix_disps = fix_disps, fix_alphas = fix_alphas,
        same_disps = same_disps, same_alphas = same_alphas)
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll2(
        data, new_params,
        weights_and_nodes, family = "cmp",
        fix_disps = fix_disps, fix_alphas = fix_alphas,
        same_disps = same_disps, same_alphas = same_alphas)
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



# truncated quadrature rule
quad_rule <- function(n_nodes, thres = Inf, prob = 0) {
  weights_and_nodes <- gaussHermiteData(n_nodes)
  # rescale because we are approximating integral for normal density
  weights_and_nodes$x <- weights_and_nodes$x * sqrt(2)
  weights_and_nodes$w <- weights_and_nodes$w / sqrt(pi)
  trunc <- abs(weights_and_nodes$x) < thres & 
    weights_and_nodes$w > prob
  weights_and_nodes$x <- weights_and_nodes$x[trunc]
  weights_and_nodes$w <- weights_and_nodes$w[trunc] 
  return(weights_and_nodes)
}


