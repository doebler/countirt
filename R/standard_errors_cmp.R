# gradients for full 2pcmp model standard errors -------------------------------------
grad_for_se <- function(y, item_params, weights_and_nodes, data, item_offset) {
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  post_probs <- newem_estep2(data = data, 
                             item_params = y, 
                             weights_and_nodes = weights_and_nodes,
                             item_offset = item_offset)
  g <- grad_cmp_newem2(item_params = item_params,
                       PPs = post_probs,
                       weights_and_nodes = weights_and_nodes,
                       data = data,
                       item_offset = item_offset)
  return(g)
}

wrap_grad_cmp_newem2 <- function(y, PPs, weights_and_nodes, data, item_offset) {
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  grad <- grad_cmp_newem2(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    item_offset = item_offset
  )
  return(grad)
}

# gradients for 2pcmp with constant alphas ----------------------------------------
grad_for_se_same_alpha <- function(y, item_params, weights_and_nodes, data, item_offset) {
  # item params and y here only have one alpha at the start and 
  # then the remaining item parameters
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  # prep the parameters for the e-step
  alpha <- y[grepl("alpha", names(y))]
  deltas <- y[grepl("delta", names(y))]
  n_items <- length(deltas)
  log_disps <- y[grepl("log_disp", names(y))]
  item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps)
  names(item_params_samealph) <- c(
    paste0("alpha", 1:n_items),
    names(y[grepl("delta", names(y))]),
    names(y[grepl("log_disp", names(y))])
  )
  
  post_probs <- newem_estep2(data = data, 
                             item_params = item_params_samealph, 
                             weights_and_nodes = weights_and_nodes,
                             item_offset = item_offset)
  g <- grad_cmp_samealphas_newem(item_params = item_params,
                                 PPs = post_probs,
                                 weights_and_nodes = weights_and_nodes,
                                 data = data,
                                 item_offset = item_offset)
  return(g)
}

wrap_grad_cmp_samealphas_newem <- function(y, PPs, weights_and_nodes, data, item_offset) {
  # y just have one alpha and then the remaining item parameters
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  grad <- grad_cmp_samealphas_newem(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    item_offset = item_offset
  )
  return(grad)
}

# gradients for 2pcmp model with constant dispersions standard errors ------------
grad_for_se_same_disp <- function(y, item_params, weights_and_nodes, data, item_offset) {
  # y and item parameters have alphas and deltas for each item but just one
  # disp parameter
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  # prep the parameters for the e-step
  alphas <- y[grepl("alpha", names(y))]
  deltas <- y[grepl("delta", names(y))]
  n_items <- length(deltas)
  log_disp <- y[grepl("log_disp", names(y))]
  item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items))
  names(item_params_samedisp) <- c(
    names(y[grepl("alpha", names(y))]),
    names(y[grepl("delta", names(y))]),
    paste0("log_disp", 1:n_items)
  )
  
  post_probs <- newem_estep2(data = data, 
                             item_params = item_params_samedisp, 
                             weights_and_nodes = weights_and_nodes,
                             item_offset = item_offset)
  g <- grad_cmp_samedisps_newem(item_params = item_params,
                                PPs = post_probs,
                                weights_and_nodes = weights_and_nodes,
                                data = data,
                                item_offset = item_offset)
  return(g)
}

wrap_grad_cmp_samedisps_newem <- function(y, PPs, weights_and_nodes, data, item_offset) {
  # y and item parameters have alphas and deltas for each item but just one
  # disp parameter
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  grad <- grad_cmp_samedisps_newem(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    item_offset
  )
  return(grad)
}

# gradients for 2pcmp model with fixed dispersions standard errors -----------------
grad_for_se_fix_disps <- function(y, item_params, weights_and_nodes, data, fix_disps, item_offset) {
  # y and item parameters only have alphas and deltas, dispersions
  # are fixed an contained in fix_disps
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  # prep the parameters for the e-step
  item_params_fixdisps <- c(item_params, log(fix_disps))
  names(item_params_fixdisps) <- c(names(item_params), 
                                   paste0("log_disp", 1:length(fix_disps)))
  post_probs <- newem_estep2(data, item_params_fixdisps, weights_and_nodes, item_offset = item_offset)
  
  g <- grad_cmp_fixdisps_newem(item_params = item_params,
                               PPs = post_probs,
                               weights_and_nodes = weights_and_nodes,
                               data = data,
                               fix_disps = fix_disps,
                               item_offset = item_offset)
  return(g)
}

wrap_grad_cmp_fixdisps_newem <- function(y, PPs, weights_and_nodes, data, fix_disps, item_offset) {
  # y and item parameters only have alphas and deltas, dispersions
  # are fixed an contained in fix_disps
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  grad <- grad_cmp_fixdisps_newem(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    fix_disps = fix_disps,
    item_offset = item_offset
  )
  return(grad)
}

# gradients for 2pcmp model with fixed alphas standard errors -------------------
grad_for_se_fix_alphas <- function(y, item_params, weights_and_nodes, data, fix_alphas, item_offset) {
  # y and item parameters only have deltas and disps, slopes
  # are fixed and contained in fix_alphas
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  # prep the parameters for the e-step
  item_params_fixalphas <- c(fix_alphas, item_params)
  names(item_params_fixalphas) <- c(paste0("alpha", 1:length(fix_alphas)), 
                                    names(item_params))
  post_probs <- newem_estep2(data, item_params_fixalphas, weights_and_nodes, item_offset = item_offset)
  
  g <- grad_cmp_fixalphas_newem(item_params = item_params,
                                PPs = post_probs,
                                weights_and_nodes = weights_and_nodes,
                                data = data,
                                fix_alphas = fix_alphas,
                                item_offset = item_offset)
  return(g)
}

wrap_grad_cmp_fixalphas_newem <- function(y, PPs, weights_and_nodes, data, fix_alphas, item_offset) {
  # y and item parameters only have deltas and disps, slopes
  # are fixed and contained in fix_alphas
  
  # handle item_offset prep in compute_vcov, so expect non-empty vector of length n_items here
  
  grad <- grad_cmp_fixalphas_newem(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    fix_alphas = fix_alphas,
    item_offset = item_offset
  )
  return(grad)
}

# compute_vcov -------------------------------------------------------------------
compute_vcov <- function(item_params, weights_and_nodes, data,
                         same_alphas = FALSE, same_disps = FALSE,
                         fix_alphas = NULL, fix_disps = NULL,
                         item_offset = NULL) {
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, ncol(data))
  } 
  
  # computes vcov matrix with Oake's identity approximation (Chalmers, 2012)#
  
  if (is.null(fix_disps) & is.null(fix_alphas)) {
    if (!same_disps & !same_alphas) {
      # standard errors for full model
      
      # compute derivative of gradient with respect to new item params
      post_probs <- newem_estep2(data, item_params, weights_and_nodes, item_offset = item_offset)
      
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_newem2,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset
      )
    } else if (same_alphas & !same_disps) {
      # standard errors for constant alphas across items
      
      # compute derivative of gradient with respect to new item params
      
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
      
      post_probs <- newem_estep2(data, item_params_samealph, weights_and_nodes, item_offset = item_offset)
      
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_samealphas_newem,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_same_alpha,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset
      )
    } else if (!same_alphas & same_disps) {
      # standard errors for constant dispersions across items
      
      # compute derivative of gradient with respect to new item params
      
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
      
      post_probs <- newem_estep2(data, item_params_samedisp, weights_and_nodes, item_offset = item_offset)
      
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_samedisps_newem,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_same_disp,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset
      )
    }
  } else {
    # we either have fixed dispersions or fixed alphas
    if (!is.null(fix_disps)) {
      # we only have alphas and deltas, but we fix dispersions to the values provided
      # in fix_disps
      
      # prep for e step
      item_params_fixdisps <- c(item_params, log(fix_disps))
      names(item_params_fixdisps) <- c(names(item_params), 
                                       paste0("log_disp", 1:length(fix_disps)))
      
      post_probs <- newem_estep2(data, item_params_fixdisps, weights_and_nodes, item_offset = item_offset)
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_fixdisps_newem,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_disps = fix_disps,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_fix_disps,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_disps = fix_disps,
        item_offset = item_offset
      )
    } else if (!is.null(fix_alphas)) {
      # we only have deltas and disps, but we fix slopes to the values provided
      # in fix_alphas
      
      # prep for e step
      item_params_fixalphas <- c(fix_alphas, item_params)
      names(item_params_fixalphas) <- c(paste0("alpha", 1:length(fix_alphas)), 
                                        names(item_params))
      
      post_probs <- newem_estep2(data, item_params_fixalphas, weights_and_nodes, item_offset = item_offset)
      
      x <- numDeriv::jacobian(
        wrap_grad_cmp_fixalphas_newem,
        item_params,
        PPs = post_probs, 
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_alphas = fix_alphas,
        item_offset = item_offset
      )
      
      x2 <- numDeriv::jacobian(
        grad_for_se_fix_alphas,
        item_params, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        data = data,
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

se_from_vcov <- function(vcov_matrix) {
  if (any(diag(vcov_matrix) < 0)) {
    warning("Some diagonal elements of the variance-covariance matrix were negative. NAs were returned.")
  }
  return(sqrt(diag(vcov_matrix)))
}