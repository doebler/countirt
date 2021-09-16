library(rootSolve)

grad_for_se <- function(y, item_params, weights_and_nodes, data) {
  post_probs <- newem_estep2(data = data, 
                             item_params = y, 
                             weights_and_nodes = weights_and_nodes)
  g <- grad_cmp_newem2(item_params = item_params,
                       PPs = post_probs,
                       weights_and_nodes = weights_and_nodes,
                       data = data)
  return(g)
}

wrap_grad_cmp_newem2 <- function(y, PPs, weights_and_nodes, data) {
  grad <- grad_cmp_newem2(
    item_params = y,
    PPs = PPs,
    weights_and_nodes = weights_and_nodes,
    data = data
  )
  return(grad)
}

compute_vcov <- function(item_params, weights_and_nodes, data) {
  # TODO implement all the constraints
  
  # computes vcov matrix with Oake's identity approximation (Chalmers, 2012)#
  
  # compute derivative of gradient with respect to new item params
  post_probs <- newem_estep2(data, item_params, weights_and_nodes)

  
  x <- numDeriv::jacobian(
      wrap_grad_cmp_newem2,
      item_params,
      PPs = post_probs, 
      weights_and_nodes = weights_and_nodes,
      data = data
    )
    
  x2 <- numDeriv::jacobian(
      grad_for_se,
      item_params, 
      item_params = item_params,
      weights_and_nodes = weights_and_nodes,
      data = data
    )
    
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
  





