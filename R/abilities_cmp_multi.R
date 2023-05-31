# get_ability_params_multi_pp ------------------------------------------------------------
# FIXME will so far only work with GH quadrature
get_ability_params_multi_pp <- function(fit, n_nodes = NULL, truncate_grid = NULL) {
  
  # FIXME i just put the n_nodes argument and this in here because old mcirt fits
  # won't have the n_nodes and the truncate grid in the model list
  if (is.null(n_nodes)) {
    n_nodes <- fit$model$n_nodes
  }
  
  if (is.null(truncate_grid)) {
    truncate_grid <- fit$model$truncate_grid
  }
  
  # get weights and nodes
  # get nodes and weights for multivariate GH quadrature
  if (is.null(fit$model$fcov_prior)) {
    weights_and_nodes <- init.quad(
      Q = fit$model$nfactors, 
      ip = n_nodes, 
      prune = truncate_grid)
  } else {
    weights_and_nodes <- init.quad(
      Q = fit$model$nfactors, 
      prior = fit$model$fcov_prior,
      ip = n_nodes, 
      prune = truncate_grid)
  }
  # weights W are on log scale
  
  # compute post probs of the thetas with the final version of the item parameters
  post_probs <- e_step_multi(
    data = fit$model$data,
    item_params = fit$fit$params,
    n_traits = fit$model$nfactors,
    em_type = fit$model$em_type,
    weights_and_nodes = weights_and_nodes
  ) # output: a matrix with N rows (= no. of persons) and K columns (= no. of node combinations)
  
  # poit ability estimates
  thetas <- data.frame(
    theta1 = rep(NA, nrow(fit$model$data))
  )
  thetas$theta1 <- apply(post_probs, 1, function(y) {sum(y*weights_and_nodes$X[,1])})
  for (i in 2:fit$model$nfactors) {
    thetas[[paste0("theta",i)]] <- apply(post_probs, 1, function(y) {sum(y*weights_and_nodes$X[,i])})
  } # each theta is of length N
  
  # standard error estimates
  se_thetas <- data.frame(
    se_theta1 = rep(NA, nrow(fit$model$data))
  )
  se_theta1 <- outer(weights_and_nodes$X[,1], thetas$theta1, "-") # K times N
  se_theta1 <- t(se_theta1)^2 * post_probs 
  # transpose se_theta1 and then multiply element wise with posterior probs
  se_theta1 <- sqrt(apply(se_theta1, 1, sum))
  se_thetas$se_theta1 <- se_theta1
  for (i in 2:fit$model$nfactors) {
    se <- outer(weights_and_nodes$X[,i], thetas[,i], "-") # K times N
    se <- t(se)^2 * post_probs 
    # transpose se_theta1 and then multiply element wise with posterior probs
    se <- sqrt(apply(se, 1, sum))
    se_thetas[[paste0("se_theta",i)]] <- se
  }
  
  # CIs
  # FIXME i dont know yet that this work correctly
  # index_lower <- apply(post_probs, 1, function(x) {
  #   ind <- which(cumsum(x) <= .025)
  #   out <- ifelse(length(ind) == 0, 1, max(ind))
  #   return(out)
  # })
  # index_upper <- apply(post_probs, 1, function(x) {min(which(cumsum(x) >= .975))})
  # ci_thetas <- data.frame(
  #   CI_lower_theta1 = rep(NA, nrow(fit$model$data)),
  #   CI_upper_theta1 = rep(NA, nrow(fit$model$data))
  # )
  # for (i in 1:fit$model$nfactors) {
  #   ci_thetas[[paste0("CI_lower_theta",i)]] <- weights_and_nodes$X[index_lower,i]
  #   ci_thetas[[paste0("CI_upper_theta",i)]] <- weights_and_nodes$X[index_upper,i]  
  # }
  out <- cbind(thetas, se_thetas) # ci_thetas
  return(out)
}