# gradients for full 2pcmp model standard errors -------------------------------------
grad_for_se_cmp_with_cov <- function(y, item_params, weights_and_nodes, data,
                                     p_covariates, i_covariates, item_offset,
                                     i_cov_on = c("alpha", "delta", "log_disp"),
                                     which_i_cov = list(alpha="all", delta="all", log_disp="all"),
                                     p_cov_cat = TRUE,
                                     resp_patterns_matrix = NULL) {
  post_probs <- estep_cmp_with_cov(
    data = data, 
    item_params = y, 
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    which_i_cov = which_i_cov,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
    )
  g <- grad_cmp_with_cov(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    which_i_cov = which_i_cov,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  return(g)
}

wrap_grad_cmp_with_cov <- function(y, PPs, weights_and_nodes, data,
                                   p_covariates, i_covariates, item_offset,
                                   i_cov_on = c("alpha", "delta", "log_disp"),
                                   which_i_cov = list(alpha="all", delta="all", log_disp="all"),
                                   p_cov_cat = TRUE,
                                   resp_patterns_matrix = NULL) {
  grad <- grad_cmp_with_cov(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    which_i_cov = which_i_cov,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  return(grad)
}

# TODO when i've implemented which_i_cov for the constrained gradients
# also add that argument here
# gradients for 2pcmp with constant alphas ----------------------------------------
grad_for_se_cmp_samealpha_with_cov <- function(y, item_params, weights_and_nodes, data,
                                               p_covariates, i_covariates, item_offset,
                                               i_cov_on = c("delta", "log_disp"),
                                               p_cov_cat = TRUE, 
                                               resp_patterns_matrix = NULL) {
  # item params and y here only have one alpha at the start and 
  # then the remaining item parameters
  
  # prep the parameters for the e-step
  alpha <- item_params[grepl("alpha", names(item_params)) &
                         !grepl("beta", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  n_items <- ncol(data)
  log_disps <- item_params[grepl("log_disp", names(item_params)) &
                             !grepl("beta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  # if we have the constraint of same alphas, we can either have covariates
  # on just delta or just log_disp, in which case we just have betas_i,
  # or we could have covariates on delta and log_disp together (or person covariates)
  # with the constraint, we can't have item covariates on all three parameters
  if (length(i_cov_on) == 2) {
    # must be covariates on delta and log_disp together
    betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
    betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
    item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps,  
                              betas_i_delta, betas_i_logdisp)
    names(item_params_samealph) <- c(
      paste0("alpha", 1:n_items),
      names(item_params[grepl("delta", names(item_params)) & 
                          !grepl("beta", names(item_params))]),
      names(item_params[grepl("log_disp", names(item_params)) & 
                          !grepl("beta", names(item_params))]),
      names(item_params[grepl("beta_i_delta", names(item_params))]),
      names(item_params[grepl("beta_i_log_disp", names(item_params))])
    )
  } else {
    # either we have person covariates or just item covariates on one parameter, so just beta_i
    item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps, betas_p, betas_i)
    names(item_params_samealph) <- c(
      paste0("alpha", 1:n_items),
      names(item_params[grepl("delta", names(item_params))]),
      names(item_params[grepl("log_disp", names(item_params))]),
      names(item_params[grepl("beta_p", names(item_params))]),
      names(item_params[grepl("beta_i", names(item_params))])
    )
  }
  
  post_probs <- estep_cmp_with_cov(
    data = data, 
    item_params = item_params_samealph, 
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates, 
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  g <- grad_cmp_with_cov_samealphas(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates, 
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  return(g)
}

wrap_grad_cmp_samealpha_with_cov <- function(y, PPs, weights_and_nodes, data,
                                             p_covariates, i_covariates, item_offset,
                                             i_cov_on = c("delta", "log_disp"),
                                             p_cov_cat = TRUE, 
                                             resp_patterns_matrix = NULL) {
  # y just have one alpha and then the remaining item parameters
  grad <- grad_cmp_with_cov_samealphas(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates, 
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  return(grad)
}

# gradients for 2pcmp model with constant dispersions standard errors ------------
grad_for_se_cmp_samedisp_with_cov <- function(y, item_params, weights_and_nodes, data,
                                              p_covariates, i_covariates, item_offset,
                                              i_cov_on = c("alpha", "delta"),
                                              p_cov_cat = TRUE,
                                              resp_patterns_matrix = NULL) {
  # y and item parameters have alphas and deltas for each item but just one
  # disp parameter
  
  # prep the parameters for the e-step
  alphas <- item_params[grepl("alpha", names(item_params)) &
                          !grepl("beta", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  n_items <- ncol(data)
  log_disp <- item_params[grepl("log_disp", names(item_params)) &
                            !grepl("beta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  # if we have the constraint of same disps, we can either have item covariates 
  # on just alpha or just delta or both together (can't have covariates on all three)
  # (or person covariates)
  # but with the constraint, we can't have item covariates on all three parameters
  if (length(i_cov_on) == 2) {
    betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
    betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
    item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items), 
                              betas_i_alpha, betas_i_delta)
    names(item_params_samedisp) <- c(
      names(item_params[grepl("alpha", names(item_params))]),
      names(item_params[grepl("delta", names(item_params))]),
      paste0("log_disp", 1:n_items),
      names(item_params[grepl("beta_i_alpha", names(item_params))]),
      names(item_params[grepl("beta_i_delta", names(item_params))])
    )
  } else {
    # either we have person covariates or just item covariates on one parameter, so just beta_i
    item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items), betas_p, betas_i)
    names(item_params_samedisp) <- c(
      names(item_params[grepl("alpha", names(item_params))]),
      names(item_params[grepl("delta", names(item_params))]),
      paste0("log_disp", 1:n_items),
      names(item_params[grepl("beta_p", names(item_params))]),
      names(item_params[grepl("beta_i", names(item_params))])
    )
  }
  
  post_probs <- estep_cmp_with_cov(
    data = data, 
    item_params = item_params_samedisp, 
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates, 
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  g <- grad_cmp_with_cov_samedisps(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates, 
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  return(g)
}

wrap_grad_cmp_samedisp_with_cov <- function(y, PPs, weights_and_nodes, data,
                                            p_covariates, i_covariates, item_offset,
                                            i_cov_on = c("alpha", "delta"),
                                            p_cov_cat = TRUE,
                                            resp_patterns_matrix = NULL) {
  # y and item parameters have alphas and deltas for each item but just one
  # disp parameter
  
  grad <- grad_cmp_with_cov_samedisps(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates, 
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  return(grad)
}

# gradients for 2pcmp model with fixed dispersions standard errors -----------------
grad_for_se_cmp_fixdisps_with_cov <- function(y, item_params, weights_and_nodes, 
                                              data, fix_disps, item_offset,
                                              p_covariates, i_covariates,
                                              i_cov_on = c("alpha", "delta"),
                                              p_cov_cat = TRUE,
                                              resp_patterns_matrix = NULL) {
  # y and item parameters only have alphas and deltas, dispersions
  # are fixed an contained in fix_disps
  
  # prep the parameters for the e-step
  # prep for e step
  alphas <- item_params[grepl("alpha", names(item_params)) &
                          !grepl("beta", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  n_items <- ncol(data)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  # if we have the constraint of fixed disps, we can either have person covariates,
  # or item covariates on alpha or delta (then we just have betas_i) or we have
  # item covariates on alpha and delta together
  # but with the constraint, we can't have item covriates on all three parameters
  if (length(i_cov_on) == 2) {
    betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
    betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
    item_params_samedisp <- c(alphas, deltas, log(fix_disps), 
                              betas_i_alpha, betas_i_delta)
    names(item_params_samedisp) <- c(
      names(item_params[grepl("alpha", names(item_params))]),
      names(item_params[grepl("delta", names(item_params))]),
      paste0("log_disp", 1:n_items),
      names(item_params[grepl("beta_i_alpha", names(item_params))]),
      names(item_params[grepl("beta_i_delta", names(item_params))])
    )
  } else {
    # either we have person covariates or just item covariates on one parameter, so just beta_i
    item_params_samedisp <- c(alphas, deltas, log(fix_disps), betas_p, betas_i)
    names(item_params_samedisp) <- c(
      names(item_params[grepl("alpha", names(item_params))]),
      names(item_params[grepl("delta", names(item_params))]),
      paste0("log_disp", 1:n_items),
      names(item_params[grepl("beta_p", names(item_params))]),
      names(item_params[grepl("beta_i", names(item_params))])
    )
  }
  
  post_probs <- estep_cmp_with_cov(
    data = data, 
    item_params = item_params_fixdisps, 
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  
  g <- grad_cmp_with_cov_fixdisps(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    fix_disps = fix_disps,
    item_offset = item_offset
    )
  return(g)
}

wrap_grad_cmp_fixdisps_with_cov <- function(y, PPs, weights_and_nodes, 
                                            data, fix_disps, item_offset,
                                            p_covariates, i_covariates,
                                            i_cov_on = c("alpha", "delta"),
                                            p_cov_cat = TRUE,
                                            resp_patterns_matrix = NULL) {
  # y and item parameters only have alphas and deltas, dispersions
  # are fixed an contained in fix_disps
  
  grad <- grad_cmp_with_cov_fixdisps(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    fix_disps = fix_disps,
    item_offset = item_offset
  )
  return(grad)
}

# gradients for 2pcmp model with fixed alphas standard errors -------------------
grad_for_se_cmp_fixalphas_with_cov <- function(y, item_params, weights_and_nodes, 
                                   data, fix_alphas, item_offset,
                                   p_covariates, i_covariates,
                                   i_cov_on = c("delta", "log_disp"),
                                   p_cov_cat = TRUE, 
                                   resp_patterns_matrix = NULL) {
  # y and item parameters only have deltas and disps, slopes
  # are fixed and contained in fix_alphas
  
  # prep for e step
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  n_items <- ncol(data)
  log_disps <- item_params[grepl("log_disp", names(item_params)) &
                             !grepl("beta", names(item_params))]
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  # if we have the constraint of fixed alphas, we can either have covariates
  # on just delta or just log_disp, in which case we just have betas_i,
  # or we could have covariates on delta and log_disp together (or person covariates)
  # but with the constraint, we can't have item covaraites on all three parameters
  if (length(i_cov_on) == 2) {
    # must be covariates on delta and log_disp together
    betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
    betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
    item_params_samealph <- c(fix_alphas, deltas, log_disps,  
                              betas_i_delta, betas_i_logdisp)
    names(item_params_samealph) <- c(
      paste0("alpha", 1:n_items),
      names(item_params[grepl("delta", names(item_params)) & 
                          !grepl("beta", names(item_params))]),
      names(item_params[grepl("log_disp", names(item_params))]),
      names(item_params[grepl("beta_i_delta", names(item_params))]),
      names(item_params[grepl("beta_i_log_disp", names(item_params))])
    )
  } else {
    # either we have person covariates or just item covariates on one parameter, so just beta_i
    item_params_samealph <- c(fix_alphas, deltas, log_disps, betas_p, betas_i)
    names(item_params_samealph) <- c(
      paste0("alpha", 1:n_items),
      names(item_params[grepl("delta", names(item_params))]),
      names(item_params[grepl("log_disp", names(item_params))]),
      names(item_params[grepl("beta_p", names(item_params))]),
      names(item_params[grepl("beta_i", names(item_params))])
    )
  }
  
  post_probs <- estep_cmp_with_cov(
    data = data, 
    item_params = item_params_fixalphas,
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    item_offset = item_offset
  )
  
  g <- grad_cmp_with_cov_fixalphas(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    fix_alphas = fix_alphas,
    item_offset = item_offset
  )
  return(g)
}

wrap_grad_cmp_fixalphas_with_cov <- function(y, PPs, weights_and_nodes, 
                                             data, fix_alphas, item_offset,
                                             p_covariates, i_covariates,
                                             i_cov_on = c("delta", "log_disp"),
                                             p_cov_cat = TRUE, 
                                             resp_patterns_matrix = NULL) {
  # y and item parameters only have deltas and disps, slopes
  # are fixed and contained in fix_alphas
  
  grad <- grad_cmp_with_cov_fixalphas(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix,
    fix_alphas = fix_alphas,
    item_offset = item_offset
  )
  return(grad)
}

# compute_vcov_with_cov -------------------------------------------------------------------
compute_vcov_with_cov <- function(item_params, weights_and_nodes, data,
                         p_covariates, i_covariates,
                         i_cov_on = c("alpha", "delta", "log_disp"),
                         which_i_cov = list(alpha="all", delta="all", log_disp="all"),
                         p_cov_cat = TRUE,
                         resp_patterns_matrix = NULL,
                         same_alphas = FALSE, same_disps = FALSE,
                         fix_alphas = NULL, fix_disps = NULL,
                         item_offset = NULL) {
  
  # computes vcov matrix with Oake's identity approximation (Chalmers, 2012)#
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, ncol(data))
  } 
  
  if (is.null(fix_disps) & is.null(fix_alphas)) {
    if (!same_disps & !same_alphas) {
      # standard errors for full model
      
      # compute derivative of gradient with respect to new item params
      post_probs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params, 
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        which_i_cov = which_i_cov,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_with_cov,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        which_i_cov = which_i_cov,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_cmp_with_cov,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        which_i_cov = which_i_cov,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
    } else if (same_alphas & !same_disps) {
      # standard errors for constant alphas across items
      
      # TODO add which_i_cov argument once i've implemented that above and
      # in the gradients for the constrained versions
      
      # compute derivative of gradient with respect to new item params
      
      # prep for e step
      alpha <- item_params[grepl("alpha", names(item_params)) &
                             !grepl("beta", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params)) &
                              !grepl("beta", names(item_params))]
      n_items <- ncol(data)
      log_disps <- item_params[grepl("log_disp", names(item_params)) &
                                 !grepl("beta", names(item_params))]
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      
      # if we have the constraint of same alphas, we can either have covariates
      # on just delta or just log_disp, in which case we just have betas_i,
      # or we could have covariates on delta and log_disp together (or person covariates)
      # with the constraint, we can't have item covariates on all three parameters
      if (length(i_cov_on) == 2) {
        # must be covariates on delta and log_disp together
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
        item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps,  
                                  betas_i_delta, betas_i_logdisp)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params)) & 
                              !grepl("beta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params)) & 
                              !grepl("beta", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))]),
          names(item_params[grepl("beta_i_log_disp", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps, betas_p, betas_i)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params))]),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }
      
      post_probs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_samealph, 
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_samealpha_with_cov,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_cmp_samealpha_with_cov,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
    } else if (!same_alphas & same_disps) {
      # standard errors for constant dispersions across items
      
      # compute derivative of gradient with respect to new item params
      
      # prep the parameters for the e-step
      # prep the parameters for the e-step
      alphas <- item_params[grepl("alpha", names(item_params)) &
                              !grepl("beta", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params)) &
                              !grepl("beta", names(item_params))]
      n_items <- ncol(data)
      log_disp <- item_params[grepl("log_disp", names(item_params)) &
                                !grepl("beta", names(item_params))]
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      
      # if we have the constraint of same disps, we can either have item covariates 
      # on just alpha or just delta or both together (can't have covariates on all three)
      # (or person covariates)
      # but with the constraint, we can't have item covariates on all three parameters
      if (length(i_cov_on) == 2) {
        betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items), 
                                  betas_i_alpha, betas_i_delta)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_i_alpha", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items), betas_p, betas_i)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }
      
      post_probs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_samedisp, 
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_samedisp_with_cov,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_cmp_samedisp_with_cov,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
    }
  } else {
    # we either have fixed dispersions or fixed alphas
    if (!is.null(fix_disps)) {
      # we only have alphas and deltas, but we fix dispersions to the values provided
      # in fix_disps
      
      # prep for e step
      alphas <- item_params[grepl("alpha", names(item_params)) &
                              !grepl("beta", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params)) &
                              !grepl("beta", names(item_params))]
      n_items <- ncol(data)
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      
      # if we have the constraint of fixed disps, we can either have person covariates,
      # or item covariates on alpha or delta (then we just have betas_i) or we have
      # item covariates on alpha and delta together
      # but with the constraint, we can't have item covriates on all three parameters
      if (length(i_cov_on) == 2) {
        betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        item_params_samedisp <- c(alphas, deltas, log(fix_disps), 
                                  betas_i_alpha, betas_i_delta)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_i_alpha", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samedisp <- c(alphas, deltas, log(fix_disps), betas_p, betas_i)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }
      
      post_probs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_fixdisps, 
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_fixdisps_with_cov,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_disps = fix_disps,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_cmp_fixdisps_with_cov,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_disps = fix_disps,
        item_offset = item_offset
      )
    } else if (!is.null(fix_alphas)) {
      # we only have deltas and disps, but we fix slopes to the values provided
      # in fix_alphas
      
      # prep for e step
      deltas <- item_params[grepl("delta", names(item_params)) &
                              !grepl("beta", names(item_params))]
      n_items <- ncol(data)
      log_disps <- item_params[grepl("log_disp", names(item_params)) &
                                 !grepl("beta", names(item_params))]
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      
      # if we have the constraint of fixed alphas, we can either have covariates
      # on just delta or just log_disp, in which case we just have betas_i,
      # or we could have covariates on delta and log_disp together (or person covariates)
      # but with the constraint, we can't have item covaraites on all three parameters
      if (length(i_cov_on) == 2) {
        # must be covariates on delta and log_disp together
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
        item_params_samealph <- c(fix_alphas, deltas, log_disps,  
                                  betas_i_delta, betas_i_logdisp)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params)) & 
                              !grepl("beta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))]),
          names(item_params[grepl("beta_i_log_disp", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samealph <- c(fix_alphas, deltas, log_disps, betas_p, betas_i)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params))]),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }
      
      post_probs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_fixalphas, 
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        item_offset = item_offset
      )
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_fixalphas_with_cov,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_alphas = fix_alphas,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_cmp_fixalphas_with_cov,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_alphas = fix_alphas,
        item_offset = item_offset
      )
    }
  }
  
  # compute Hessian with Oake's identity
  hess_matrix <- x + x2
  
  # ensure Hessian is resulting is symmetrical
  below_diag <- hess_matrix[lower.tri(hess_matrix)]
  above_diag <- t(hess_matrix)[lower.tri(hess_matrix)]
  avg_off_diag <- (below_diag + above_diag) / 2
  hess_matrix_symm <- diag(diag(hess_matrix))
  hess_matrix_symm[lower.tri(hess_matrix_symm)] <- avg_off_diag
  hess_matrix_symm <- t(hess_matrix_symm)
  hess_matrix_symm[lower.tri(hess_matrix_symm)] <- avg_off_diag
  
  # compute variance-covariance matrix of item parameter estimators 
  # from Hessian
  vcov_matrix <- solve(-hess_matrix_symm)
  
  if (min(eigen(vcov_matrix)$values) <= 0) {
    warning("Variance-covariance matrix is not positive semi-definite.")
  }
  
  return(vcov_matrix)
}

