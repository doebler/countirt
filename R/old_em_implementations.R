# old em code -----------------------------------------------------------------

# this is old code from previous em implementations based on implementing
# em with summary statistics
# it is not tested code, and they are just internal functions, just
# to keep as an archive and if one wants to implement the em algorithm
# based on summary stats
# there might still be mistakes in the implementations, as this is just
# archived code

# -----------------------------------------------------------------------------

# grid_W_long <- as.vector(grid_W)
# grid_R_long <- as.vector(grid_R)
# grid_W_long_nona <- ifelse(is.na(grid_W_long), 
#                            max(grid_W_long, na.rm = TRUE), 
#                            grid_W_long)
# 
# grid_R_long_nona <- ifelse(is.na(grid_R_long), 
#                            max(grid_R_long, na.rm = TRUE), 
#                            grid_R_long)

# functions for e step

# @param dens: densities for response vector to one item j (of all N participants), 
#                   summing over this vector = summing over all participants
# @param quad_weights: Quadrature weight for node k
# @return: vector of posterior probabilities (for all participants) for node k and one item j
comp_post_prob_k <- function(dens, quad_weight_k, marg_prob) {
  post_prob_k <- (dens * quad_weight_k) / marg_prob
  return(post_prob_k)
}

vlogFactorial <- function(n) {
  out <- numeric(length(n))
  for (i in 1:length(n)) {
    out[i] <- logFactorial(n[i])
  }
  return(out)
}
vlogFactorial <- Vectorize(vlogFactorial)

run_e_step2 <- function(data, item_params, weights_and_nodes, family, 
                        fix_disps = NULL, max_iter = 1, interp_method = "bicubic") {
  n_items <- ncol(data)
  n_nodes <- length(weights_and_nodes$x)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  if (family == "cmp") {
    if (is.null(fix_disps)) {
      log_disps <- item_params[grepl("log_disp", names(item_params))]
      disps <- exp(log_disps)
    } else {
      disps <- fix_disps
    }
  } else if (family == "poisson") {
    # so that a function below can work for poisson case
    disps <- rep(NA, length(n_items)) 
  }
  
  if (family == "poisson") {
    e_values <- vector(mode = "list", length = n_items)
    for(j in 1:n_items) {
      lambda_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j]) 
      marg_prob <- numeric(nrow(data))
      K <- length(lambda_j)
      N <- nrow(data)
      lambda_j_matrix <- matrix(rep(lambda_j, each = N), 
                                nrow = N, ncol = K, byrow = FALSE)
      # now when we use apply columnwise we have as many lambdas from the matrix as
      # we have responses to item j, so we can feed the response vector into the dcmp
      # function and only need to rep dispersions for length of response vector to item j
      dens_i_eachnode <- apply(lambda_j_matrix, 2, function(x){dpois(data[ , j], x)})
      # in dens_i_eachnode columns represent one node each; (rows = no. of responses)
      # so here we want to get a person- (but not node-)specific marginal likelihood
      # we therefore need to go rowwise so we still have one value for each person
      # and need to sum across all nodes = all columns
      marg_prob <- apply(dens_i_eachnode, 1, function(x){sum(x*weights_and_nodes$w)})
      # compute posterior probabilties for all nodes and persons (for item j) in
      # the shape of a matrix with dimensions persons (rows) + nodes (columns)
      post_probs <- computepp_allnodes(dens_i_eachnode, weights_and_nodes$w, marg_prob)
      r_j <- computer_allnodes(data[ , j], post_probs)
      f_j <- computef_allnodes(post_probs)
      h_j <- NULL
      e_values[[j]] <- list(r_j = r_j, f_j = f_j, h_j = h_j)
    }
    r <- list.cbind(lapply(e_values, function(x){x$r_j}))
    f <- list.cbind(lapply(e_values, function(x){x$f_j}))
    h <- list.cbind(lapply(e_values, function(x){x$h_j}))
  } else if (family == "cmp") {
    if (interp_method == "bicubic") {
      suff_stats <- e_values_cpp(data = as.matrix(data),
                                 alphas = alphas, 
                                 deltas = deltas, 
                                 disps = disps, 
                                 nodes = weights_and_nodes$x,
                                 weights = weights_and_nodes$w,
                                 grid_mus = grid_mus,
                                 grid_nus = grid_nus, 
                                 grid_logZ_long = grid_logZ_long,
                                 grid_log_lambda_long = grid_log_lambda_long,
                                 max_mu = 150,
                                 min_mu = 0.001)
    } else {
      # interpolation method is linear
      suff_stats <- e_values_cpp_lininterp(data = as.matrix(data),
                                           alphas = alphas, 
                                           deltas = deltas, 
                                           disps = disps, 
                                           nodes = weights_and_nodes$x,
                                           weights = weights_and_nodes$w,
                                           grid_mus = grid_mus,
                                           grid_nus = grid_nus, 
                                           grid_logZ_long = grid_logZ_long,
                                           grid_log_lambda_long = grid_log_lambda_long, 
                                           max_mu = 200)
    }
    
    r <- suff_stats[1:n_nodes,]
    f <- suff_stats[(n_nodes+1):(2*n_nodes),]
    h <- suff_stats[(2*n_nodes+1):(3*n_nodes),]
  }
  
  out = list(r = r, f = f, h = h)
  return(out)
}

#  gradient functions ---------------------------------------------------------------------

grad_e_ll_pois <- function(item_params, e_values, weights_and_nodes) {
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(alphas)
  
  grads <- numeric(length(item_params))
  for(j in 1:length(item_params)) {
    if(grepl("alpha", names(item_params)[j])) {
      # compute lambdas for current item parameter and all nodes
      # we need the expected counts for the current item j at each of node k
      lambda_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
      # lambda_j is a vector with k elements here, as is r_j and f_j
      grads[j] <- sum(weights_and_nodes$x * (e_values[[j]][["r_j"]] -
                                               lambda_j * e_values[[j]][["f_j"]]))
    } else {
      # compute lambdas for current item parameter and all nodes
      lambda_j <- exp(alphas[j - n_items] * weights_and_nodes$x + deltas[j - n_items])
      grads[j] <- sum(e_values[[j - n_items]][["r_j"]] - lambda_j * e_values[[j - n_items]][["f_j"]])
    }
  }
  return(grads)
}

grad_e_ll_pois_fixa <- function(item_params, e_values, weights_and_nodes, alphas) {
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(deltas)
  
  grads <- numeric(length(item_params))
  for(j in 1:length(item_params)) {
    # compute lambdas for current item parameter and all nodes
    lambda_j <- exp(alphas[j]*weights_and_nodes$x + deltas[j])
    grads[j] <- sum(e_values[[j]][["r_j"]] - lambda_j * e_values[[j]][["f_j"]])
  }
  return(grads)
}

grad_cmp <- function(item_params, e_values, weights_and_nodes) {
  # print(item_params)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  grads <- grad_cmp_cpp(alphas = alphas, deltas = deltas, disps = disps, 
                        nodes = weights_and_nodes$x,
                        weights = weights_and_nodes$w,
                        r = e_values$r, f = e_values$f, h = e_values$h,
                        grid_mus = grid_mus, grid_nus = grid_nus, 
                        grid_cmp_var_long = grid_cmp_var_long,
                        grid_log_lambda_long = grid_log_lambda_long,
                        grid_logZ_long = grid_logZ_long,
                        grid_R_long = grid_R_long_nona,
                        grid_W_long = grid_W_long_nona,
                        max_mu = 200, 
                        min_mu = 0.001)
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

grad_cmp_fixdisps <- function(item_params, e_values, weights_and_nodes, fix_disps) {
  # print(item_params)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  disps <- fix_disps
  n_items <- length(alphas)
  
  grads <- grad_cmp_fixdisps_cpp(alphas = alphas, 
                                 deltas = deltas, 
                                 disps = disps, 
                                 nodes = weights_and_nodes$x,
                                 r = e_values$r, 
                                 f = e_values$f,
                                 grid_mus = grid_mus, 
                                 grid_nus = grid_nus, 
                                 grid_cmp_var_long = grid_cmp_var_long,
                                 grid_log_lambda_long = grid_log_lambda_long,
                                 max_mu = 200)
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

grad_cmp_fixalphas <- function(item_params, e_values, weights_and_nodes, fix_alphas) {
  # print(item_params)
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- length(deltas)
  
  grads <- grad_cmp_fixalphas_cpp(alphas = alphas, 
                                  deltas = deltas, 
                                  disps = disps, 
                                  nodes = weights_and_nodes$x,
                                  r = e_values$r, 
                                  f = e_values$f,
                                  h = e_values$h, 
                                  grid_mus = grid_mus, 
                                  grid_nus = grid_nus, 
                                  grid_cmp_var_long = grid_cmp_var_long,
                                  grid_log_lambda_long = grid_log_lambda_long,
                                  grid_logZ_long = grid_logZ_long,
                                  max_mu = 200)
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

grad_cmp_fixalphas_num <- function(item_params, e_values, weights_and_nodes, fix_alphas) {
  # print(item_params)
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- length(deltas)
  
  grads <- numDeriv::grad(ell_cmp_fixalphas, item_params,
                          e_values = e_values,
                          weights_and_nodes = weights_and_nodes,
                          fix_alphas = fix_alphas)
  
  #print(grads)
  return(grads)
}

grad_cmp_logdisps <- function(log_disps, e_values, weights_and_nodes, alphas_deltas) {
  # print(item_params)
  alphas <- alphas_deltas[grepl("alpha", names(alphas_deltas))]
  deltas <- alphas_deltas[grepl("delta", names(alphas_deltas))]
  disps <- exp(log_disps)
  
  grads <- grad_cmp_logdisps_cpp(alphas = alphas, 
                                 deltas = deltas, 
                                 disps = disps, 
                                 nodes = weights_and_nodes$x,
                                 weights = weights_and_nodes$w,
                                 r = e_values$r, 
                                 f = e_values$f,
                                 h = e_values$h,
                                 grid_mus = grid_mus, 
                                 grid_nus = grid_nus, 
                                 grid_cmp_var_long = grid_cmp_var_long,
                                 grid_log_lambda_long = grid_log_lambda_long,
                                 grid_logZ_long = grid_logZ_long,
                                 max_mu = 200)
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

grad_cmp_samealpha <- function(item_params, e_values, weights_and_nodes) {
  # print(item_params)
  alpha <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  grads <- grad_cmp_samealpha_cpp(alpha = alpha, 
                                  deltas = deltas, 
                                  disps = disps, 
                                  nodes = weights_and_nodes$x,
                                  weights = weights_and_nodes$w,
                                  r = e_values$r, 
                                  f = e_values$f,
                                  h = e_values$h,
                                  grid_mus = grid_mus, 
                                  grid_nus = grid_nus, 
                                  grid_cmp_var_long = grid_cmp_var_long,
                                  grid_log_lambda_long = grid_log_lambda_long,
                                  grid_logZ_long = grid_logZ_long,
                                  max_mu = 200)
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

grad_cmp_samedisp <- function(item_params, e_values, weights_and_nodes,
                              interp_method = "bicubic") {
  # print(item_params)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disp <- item_params[grepl("log_disp", names(item_params))]
  disp <- exp(log_disp)
  
  if (interp_method == "bicubic") {
    grads <- grad_cmp_samedisp_cpp(alphas = alphas, 
                                   deltas = deltas, 
                                   disp = disp, 
                                   nodes = weights_and_nodes$x,
                                   weights = weights_and_nodes$w,
                                   r = e_values$r, 
                                   f = e_values$f,
                                   h = e_values$h,
                                   grid_mus = grid_mus, 
                                   grid_nus = grid_nus, 
                                   grid_cmp_var_long = grid_cmp_var_long,
                                   grid_log_lambda_long = grid_log_lambda_long,
                                   grid_logZ_long = grid_logZ_long,
                                   grid_W_long = grid_W_long_nona,
                                   grid_R_long = grid_R_long_nona,
                                   max_mu = 150,
                                   min_mu = 0.001)
  } else {
    grads <- grad_cmp_samedisp_cpp_lininterp(alphas = alphas, 
                                             deltas = deltas, 
                                             disp = disp, 
                                             nodes = weights_and_nodes$x,
                                             weights = weights_and_nodes$w,
                                             r = e_values$r, 
                                             f = e_values$f,
                                             h = e_values$h,
                                             grid_mus = grid_mus, 
                                             grid_nus = grid_nus, 
                                             grid_cmp_var_long = grid_cmp_var_long,
                                             grid_log_lambda_long = grid_log_lambda_long,
                                             grid_W_long = grid_W_long_nona,
                                             grid_R_long = grid_R_long_nona,
                                             max_mu = 200)
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  #print(grads)
  return(grads)
}

# -----------------------------------------------------------------------------------------


expect_ll_pois <- function(item_params, e_values, weights_and_nodes) {
  # print(item_params)
  # s <<- s + 1
  # print(s)
  #e_values[[j - 2*n_items]][["h_j"]]
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  n_items <- length(alphas)
  out <- 0
  
  for (j in 1:n_items) {
    mu_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
    out <- out + sum(
      e_values[[j]][["r_j"]] * log(mu_j) -
        e_values[[j]][["f_j"]] * mu_j
    )
  }
  
  return(out)
}

ell_cmp <- function(item_params, e_values, weights_and_nodes, fix_disps = NULL) {
  # print(item_params)
  # s <<- s + 1
  # print(s)
  
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  if (is.null(fix_disps)) {
    log_disps <- item_params[grepl("log_disp", names(item_params))]
    disps <- exp(log_disps)
  } else {
    disps <- fix_disps
  }
  
  out <- ell_cmp_cpp(alphas = alphas, deltas = deltas, disps = disps, 
                     nodes = weights_and_nodes$x,
                     r = e_values$r, f = e_values$f, h = e_values$h,
                     grid_mus = grid_mus, grid_nus = grid_nus, 
                     grid_logZ_long = grid_logZ_long,
                     grid_log_lambda_long = grid_log_lambda_long,
                     max_mu = 200)
  
  return(out)
}

ell_cmp_fixalphas <- function(item_params, e_values, weights_and_nodes, fix_alphas) {
  # print(item_params)
  # s <<- s + 1
  # print(s)
  
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  out <- ell_cmp_cpp(alphas = alphas, deltas = deltas, disps = disps, 
                     nodes = weights_and_nodes$x,
                     r = e_values$r, f = e_values$f, h = e_values$h,
                     grid_mus = grid_mus, grid_nus = grid_nus, 
                     grid_logZ_long = grid_logZ_long,
                     grid_log_lambda_long = grid_log_lambda_long,
                     max_mu = 200)
  
  return(out)
}

# em cycle --------------------------------------------------------------------------------

# @param data: matrix or dataframe with as many columns as items and as many rows as
#              participants
em_cycle <- function(data, item_params, weights_and_nodes, family, 
                     fix_disps = NULL, cmp_two_step = FALSE,
                     same_disp = FALSE, same_alpha = FALSE,
                     fix_alphas = NULL,
                     ctol_maxstep = 1e-8, m_method = "NR",
                     grad_interp_method = "bicubic") {
  # ctol default wie in multiroot
  
  if (family == "poisson") {
    if (is.null(fix_alphas)) {
      # then we have alphas and deltas
      # e step
      e_values <- run_e_step(data, item_params, weights_and_nodes, "poisson")
      
      # m step 
      if (m_method == "NR") {
        new_item_params <- multiroot(
          f = grad_e_ll_pois,
          start = item_params,
          ctol = ctol_maxstep,
          e_values = e_values,
          weights_and_nodes = weights_and_nodes
        )$root
      } else if (m_method == "nleqslv") {
        new_item_params <- nleqslv(
          x = item_params,
          fn = grad_e_ll_pois,
          e_values = e_values,
          weights_and_nodes = weights_and_nodes,
          #          method = "Newton",
          control = list(xtol = ctol_maxstep)
        )$x
      }
      
    } else {
      # in this case we only have deltas, so fix alphas at one
      # e step
      item_params_fixa <- c(fix_alphas, item_params)
      names(item_params_fixa) <- c(paste0("alpha", 1:ncol(data)), names(item_params))
      e_values <- run_e_step(data, item_params_fixa, weights_and_nodes, "poisson")
      
      # m step
      new_item_params <- multiroot(
        f = grad_e_ll_pois_fixa,
        start = item_params,
        ctol = ctol_maxstep,
        e_values = e_values,
        weights_and_nodes = weights_and_nodes, 
        alphas = fix_alphas
      )$root
    }
  } else if (family == "cmp") {
    if (is.null(fix_disps) & is.null(fix_alphas)) {
      # we don't have fixed dispersions 
      # we could then want the same dispersion across items
      if (same_disp) {
        alphas <- item_params[grepl("alpha", names(item_params))]
        n_items <- length(alphas)
        deltas <- item_params[grepl("delta", names(item_params))]
        alphas_deltas <- c(alphas, deltas)
        names(alphas_deltas) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))])
        )
        log_disp <- item_params[grepl("log_disp", names(item_params))]
        item_params_samedisp <- c(alphas_deltas, rep(log_disp, n_items))
        names(item_params_samedisp) <- c(
          names(alphas_deltas),
          paste0("log_disp", 1:n_items)
        )
        
        # e-step
        e_values <- run_e_step2(data, item_params_samedisp, weights_and_nodes, 
                                "cmp", interp_method = grad_interp_method)
        
        # m-step
        new_item_params <- nleqslv(
          x = item_params,
          fn = grad_cmp_samedisp,
          e_values = e_values,
          weights_and_nodes = weights_and_nodes,
          interp_method = grad_interp_method,
          control = list(xtol = ctol_maxstep)
        )$x
        
      } else if (same_alpha) {
        # or the same discrimination across items
        alpha <- item_params[grepl("alpha", names(item_params))]
        deltas <- item_params[grepl("delta", names(item_params))]
        n_items <- length(deltas)
        log_disps <- item_params[grepl("log_disp", names(item_params))]
        item_params_samealpha <- c(rep(alpha, n_items), deltas, log_disps)
        names(item_params_samealpha) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params))])
        )
        
        # e-step
        e_values <- run_e_step2(data, item_params_samealpha, weights_and_nodes, "cmp")
        
        # m-step
        new_item_params <- nleqslv(
          x = item_params,
          fn = grad_cmp_samealpha,
          e_values = e_values,
          weights_and_nodes = weights_and_nodes,
          control = list(xtol = ctol_maxstep)
        )$x
      } else {
        # we want to estimate dispersions and discriminations for all items,
        # so we want the option between
        # "one-step" and "two-step" maximization
        # e step
        e_values <- run_e_step2(data, item_params, weights_and_nodes, "cmp")
        
        # m step
        if (cmp_two_step) {
          # first treat dispersions as fixed and optimize for item parameters
          alphas <- item_params[grepl("alpha", names(item_params))]
          deltas <- item_params[grepl("delta", names(item_params))]
          alphas_deltas <- c(alphas, deltas)
          names(alphas_deltas) <- c(
            names(item_params[grepl("alpha", names(item_params))]),
            names(item_params[grepl("delta", names(item_params))])
          )
          log_disps <- item_params[grepl("log_disp", names(item_params))]
          names(log_disps) <- names(item_params[grepl("log_disp", names(item_params))])
          disps <- exp(log_disps)
          
          new_alphas_deltas <- nleqslv(
            x = alphas_deltas,
            fn = grad_cmp_fixdisps,
            e_values = e_values,
            weights_and_nodes = weights_and_nodes,
            fix_disps = disps,
            control = list(xtol = ctol_maxstep)
          )$x
          
          # run another e-step
          e_values <- run_e_step2(data, c(new_alphas_deltas, log_disps),
                                  weights_and_nodes, "cmp")
          
          # then treat item parameters as fixed and optimize for log dispersions
          new_log_disps <- nleqslv(
            x = log_disps,
            fn = grad_cmp_logdisps,
            e_values = e_values,
            weights_and_nodes = weights_and_nodes,
            alphas_deltas = alphas_deltas,
            control = list(xtol = ctol_maxstep)
          )$x
          
          new_item_params <- c(new_alphas_deltas, new_log_disps)
          names(new_item_params) <- c(names(new_alphas_deltas), names(new_log_disps))
          
        } else {
          if (m_method == "NR") {
            new_item_params <- multiroot(
              f = grad_cmp,
              start = item_params,
              ctol = ctol_maxstep,
              e_values = e_values,
              weights_and_nodes = weights_and_nodes
            )$root
          } else if (m_method == "optim") {
            new_item_params <- optim(
              par = item_params,
              fn = ell_cmp,
              gr = grad_cmp,
              e_values = e_values,
              weights_and_nodes = weights_and_nodes,
              method = "BFGS",
              control = list(fnscale = -1, reltol = ctol_maxstep)
            )$par
          } else if (m_method == "nleqslv") {
            new_item_params <- nleqslv(
              x = item_params,
              fn = grad_cmp,
              e_values = e_values,
              weights_and_nodes = weights_and_nodes,
              #          method = "Newton",
              control = list(xtol = ctol_maxstep)
            )$x
          }
        }
      }
    } else if (!is.null(fix_disps)) {
      # dispersions are fixed, so only optimize item parameters
      # e step
      params_estep <- c(item_params, fix_disps)
      names(params_estep) <- c(names(item_params), paste0("disp", 1:length(fix_disps)))
      e_values <- run_e_step2(data, params_estep, weights_and_nodes, "cmp", fix_disps)
      
      # m step
      if (m_method == "NR") {
        new_item_params <- multiroot(
          f = grad_cmp_fixdisps,
          start = item_params,
          ctol = ctol_maxstep,
          e_values = e_values,
          weights_and_nodes = weights_and_nodes,
          fix_disps = fix_disps
        )$root
      } else if (m_method == "nleqslv") {
        new_item_params <- nleqslv(
          x = item_params,
          fn = grad_cmp_fixdisps,
          e_values = e_values,
          weights_and_nodes = weights_and_nodes,
          fix_disps = fix_disps,
          #          method = "Newton",
          control = list(xtol = ctol_maxstep)
        )$x
      }
    } else if(!is.null(fix_alphas)){
      # we want to fix discrimination
      
      # e step
      params_estep <- c(fix_alphas, item_params)
      names(params_estep) <- c(paste0("alpha", 1:length(fix_alphas)), names(item_params))
      e_values <- run_e_step2(data, params_estep, weights_and_nodes, "cmp")
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_fixalphas,
        e_values = e_values,
        weights_and_nodes = weights_and_nodes,
        fix_alphas = fix_alphas,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  }
  return(new_item_params)
}

# marginal likelihood ---------------------------------------------------------------------------------------

# marg_ll2 <- function(data, item_params, weights_and_nodes, family, fix_disps = NULL,
#                      fix_alphas = NULL, interp_method = "bicubic",
#                      same_disps = FALSE, same_alphas = FALSE) {
#   n_items <- ncol(data)
#   n_persons <- nrow(data)
#   deltas <- item_params[grepl("delta", names(item_params))]
#   if (family == "poisson") {
#     if (is.null(fix_alphas)) {
#       if (same_alphas) {
#         alpha <- item_params[grepl("alpha", names(item_params))]
#         alphas <- rep(alpha, n_items) 
#       } else {
#         alphas <- item_params[grepl("alpha", names(item_params))]
#       }
#     } else {
#       alphas <- fix_alphas
#     }
#   } else if (family == "cmp") {
#     if (is.null(fix_alphas)) {
#       if (same_alphas) {
#         alpha <- item_params[grepl("alpha", names(item_params))]
#         alphas <- rep(alpha, n_items)
#       } else {
#         alphas <- item_params[grepl("alpha", names(item_params))]
#       }
#     } else {
#       alphas <- fix_alphas
#     }
#     if (is.null(fix_disps)) {
#       if (same_disps) {
#         log_disp <- item_params[grepl("log_disp", names(item_params))]
#         disps <- rep(exp(log_disp), n_items)
#       } else {
#         log_disps <- item_params[grepl("log_disp", names(item_params))]
#         disps <- exp(log_disps)
#       }
#     } else {
#       disps <- fix_disps
#     }
#   }
#   
#   
#   if (family == "poisson") {
#     # function to compute integral with quadrature over
#     f <- function(z, data, alphas, deltas) {
#       out <- 0
#       for (j in 1:n_items) {
#         lambda <- exp(alphas[j] * z + deltas[j])
#         out <- out + (dpois(data[,j], lambda, log = TRUE))
#       }
#       return(exp(out))
#     }
#     
#     marg_prob <- numeric(n_persons)
#     for (i in 1:n_persons) {
#       marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
#                              data = data[i, , drop = FALSE], 
#                              alphas = alphas, deltas = deltas)
#     }
#     ll <- sum(log(marg_prob))
#   } else if (family == "cmp") {
#     if (interp_method == "bicubic") {
#       ll <- marg_ll_cpp(data = as.matrix(data),
#                         alphas = alphas,
#                         deltas = deltas, 
#                         disps = disps, 
#                         nodes = weights_and_nodes$x,
#                         weights = weights_and_nodes$w,
#                         grid_mus = grid_mus,  
#                         grid_nus = grid_nus, 
#                         grid_logZ_long = grid_logZ_long,
#                         grid_log_lambda_long = grid_log_lambda_long,
#                         max_mu = 150,
#                         min_mu = 0.001)
#     } else {
#       # then interpolation method is linear
#       # here we don't cap mu, so we extrapolate beyond grid values
#       ll <- marg_ll_cpp_lininterp(data = as.matrix(data),
#                                   alphas = alphas,
#                                   deltas = deltas, 
#                                   disps = disps,
#                                   nodes = weights_and_nodes$x,
#                                   weights = weights_and_nodes$w,
#                                   grid_mus = grid_mus,  
#                                   grid_nus = grid_nus, 
#                                   grid_logZ_long = grid_logZ_long,
#                                   grid_log_lambda_long = grid_log_lambda_long)
#     }
#     
#   }
#   return(ll)
# }

# EM algorithm --------------------------------------------------------------------------------------------

# for convergence check evaluate observed likelihood instead of expected likelihood

run_em <- function(data, init_params, n_nodes, family, 
                   fix_disps = NULL, cmp_two_step = FALSE,
                   same_disp = FALSE, same_alpha = FALSE,
                   fix_alphas = NULL,
                   maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8, 
                   m_method = "NR", convcrit = "marglik",
                   mll_interp_method = "bicubic",
                   grad_interp_method = "bicubic") {
  # get nodes and weights for GH quadrature
  weights_and_nodes <- gaussHermiteData(n_nodes)
  weights_and_nodes$x <- weights_and_nodes$x * sqrt(2)
  weights_and_nodes$w <- weights_and_nodes$w / sqrt(pi)
  
  new_params <- init_params
  conv <- FALSE
  iter <- 1
  
  new_ll <- 0
  marg_lls <- c()
  
  print("Start estimation...")
  
  
  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    old_params <- new_params
    new_params <- em_cycle(data, old_params, weights_and_nodes, family, 
                           fix_disps, cmp_two_step, same_disp, same_alpha,
                           fix_alphas,
                           ctol_maxstep = ctol_maxstep, m_method = m_method,
                           grad_interp_method = grad_interp_method)
    #print(new_params)
    
    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll2(
        data, new_params, 
        weights_and_nodes, family, 
        fix_disps, fix_alphas,
        mll_interp_method,
        same_disps = same_disp,
        same_alphas = same_alpha)
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll2(data, new_params, 
                          weights_and_nodes, family, 
                          fix_disps, fix_alphas,
                          mll_interp_method,
                          same_disps = same_disp,
                          same_alphas = same_alpha)
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