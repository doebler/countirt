library(rootSolve)

# gradients for full 2p poisson model with covariats standard errors -------------------------------------
grad_for_se_poisson_with_cov <- function(y, item_params, weights_and_nodes, data,
                                p_covariates, i_covariates, offset,
                                i_cov_on = c("alpha", "delta")) {
  post_probs <- estep_poisson_with_cov(
    data = data, 
    item_params = y, 
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    offset = offset
  )
  g <- grad_poisson_with_cov(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    offset = offset
  )
  return(g)
}

wrap_grad_poisson_with_cov <- function(y, PPs, weights_and_nodes, data,
                              p_covariates, i_covariates, offset,
                              i_cov_on = c("alpha", "delta")) {
  grad <- grad_poisson_with_cov(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    offset = offset
  )
  return(grad)
}

# gradients for 2p poisson with constant alphas with covariates ----------------------------------------
grad_for_se_poisson_same_alpha_with_cov <- function(y, item_params, weights_and_nodes, data,
                                                    p_covariates, i_covariates, offset,
                                                    i_cov_on = "delta") {
  # item params and y here only have one alpha at the start and 
  # then the remaining item parameters
  
  # prep the parameters for the e-step
  alpha <- item_params[grepl("alpha", names(item_params))]
  item_params_samea <- c(rep(alpha, ncol(data)), item_params[-alpha])
  names(item_params_samea) <- c(paste0("alpha", 1:ncol(data)), 
                                names(item_params[-alpha]))
  
  post_probs <- estep_poisson_with_cov(
    data = data, 
    item_params = item_params_samea, 
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    offset = offset
  )
  
  g <- grad_poisson_with_cov_samealpha(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    offset = offset
  )
  return(g)
}

wrap_grad_poisson_samealpha_with_cov <- function(y, PPs, weights_and_nodes, data,
                                        p_covariates, i_covariates, offset,
                                        i_cov_on = "delta") {
  # y just have one alpha and then the remaining item parameters
  grad <- grad_poisson_with_cov_samealpha(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    offset = offset
  )
  return(grad)
}

# gradients for 2p poisson model with fixed alphas with covariates standard errors ----------------------
grad_for_se_poisson_fix_alphas_with_cov <- function(y, item_params, weights_and_nodes, 
                                           data, p_covariates, i_covariates, fix_alphas,
                                           offset,
                                           i_cov_on = "delta") {
  # y and item parameters only have deltas and disps, slopes
  # are fixed and contained in fix_alphas
  
  # prep the parameters for the e-step
  item_params_fixalphas <- c(fix_alphas, item_params)
  names(item_params_fixalphas) <- c(paste0("alpha", 1:length(fix_alphas)), 
                                    names(item_params))
  post_probs <- estep_poisson_with_cov(
    data = data, 
    item_params = item_params_fixalphas, 
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    offset = offset
  )
  
  g <- grad_poisson_with_cov_fixalphas(
    item_params = item_params,
    PPs = post_probs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    fix_alphas = fix_alphas,
    offset = offset
    )
  return(g)
}

wrap_grad_poisson_fixalphas_with_cov <- function(y, PPs, weights_and_nodes, data, 
                                        p_covariates, i_covariates, fix_alphas, offset,
                                        i_cov_on = "delta") {
  # y and item parameters only have deltas and disps, slopes
  # are fixed and contained in fix_alphas
  
  grad <- grad_poisson_with_cov_fixalphas(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    fix_alphas = fix_alphas,
    offset = offset
  )
  return(grad)
}

# compute_vcov_poisson -------------------------------------------------------------------
compute_vcov_poisson_with_cov <- function(item_params, weights_and_nodes, data,
                                          p_covariates, i_covariates, offset,
                                          i_cov_on = c("alpha", "delta"),
                                          same_alphas = FALSE, fix_alphas = NULL) {
  
  # computes vcov matrix with Oake's identity approximation (Chalmers, 2012)#
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, ncol(data))
  }
  
  if (!is.null(fix_alphas)) {
    # fix alphas to the provided values
    # e step
    item_params_fixa <- c(fix_alphas, item_params)
    names(item_params_fixa) <- c(paste0("alpha", 1:ncol(data)), names(item_params))
    
    post_probs <- estep_poisson_with_cov(
      data = data, 
      item_params =  item_params_fixa, 
      weights_and_nodes = weights_and_nodes,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      offset = offset
    )
    
    x <- numDeriv::jacobian(
      wrap_grad_poisson_fixalphas_with_cov,
      item_params,
      PPs = post_probs, 
      weights_and_nodes = weights_and_nodes,
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      fix_alphas = fix_alphas,
      offset = offset
    )
    
    x2 <- numDeriv::jacobian(
      grad_for_se_poisson_fix_alphas_with_cov,
      item_params, 
      item_params = item_params,
      weights_and_nodes = weights_and_nodes,
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      fix_alphas = fix_alphas,
      offset = offset
    )
  } else if (same_alphas) {
    # fit the model with estimating one same alpha for all item
    # e step
    alpha <- item_params[grepl("alpha", names(item_params))]
    item_params_samea <- c(rep(alpha, ncol(data)), 
                           item_params[!grepl("alpha", names(item_params))])
    names(item_params_samea) <- c(paste0("alpha", 1:ncol(data)), 
                                  names(item_params)[!grepl("alpha", names(item_params))])
    
    post_probs <- estep_poisson_with_cov(
      data = data,
      item_params = item_params_samea,
      weights_and_nodes = weights_and_nodes,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      offset = offset
    )
    
    x <- numDeriv::jacobian(
      wrap_grad_poisson_samealpha_with_cov,
      item_params,
      PPs = post_probs, 
      weights_and_nodes = weights_and_nodes,
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      offset = offset
    )
    
    x2 <- numDeriv::jacobian(
      grad_for_se_poisson_same_alpha_with_cov,
      item_params, 
      item_params = item_params,
      weights_and_nodes = weights_and_nodes,
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      offset = offset
    )
    
  } else {
    # fit a full two parameter model
    # e step
    post_probs <- estep_poisson_with_cov(
      data = data, 
      item_params = item_params, 
      weights_and_nodes = weights_and_nodes,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      offset = offset
    )
    
    x <- numDeriv::jacobian(
      wrap_grad_poisson_with_cov,
      item_params,
      PPs = post_probs, 
      weights_and_nodes = weights_and_nodes,
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      offset = offset
    )
    
    x2 <- numDeriv::jacobian(
      grad_for_se_poisson_with_cov,
      item_params, 
      item_params = item_params,
      weights_and_nodes = weights_and_nodes,
      data = data,
      p_covariates = p_covariates,
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      offset = offset
    )
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

  





