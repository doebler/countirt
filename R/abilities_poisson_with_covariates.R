# get_ability_params_poisson_pp_with_cov ------------------------------------------------------------

get_ability_params_poisson_pp_with_cov <- function(data, item_params, n_nodes,
                                                   p_covariates, i_covariates,
                                                   i_cov_on = c("alpha", "delta"),
                                                   same_alphas = FALSE, fix_alphas = NULL,
                                                   thres = Inf, prob = 0) {
  
  # get weights and nodes
  weights_and_nodes <- quad_rule(n_nodes, thres = thres,prob = prob)
  
  # if we have any constraints, adjust item parameters accordingly
  if (same_alphas) {
    alpha <- item_params[grepl("alpha", names(item_params))]
    item_params_estep <- c(rep(alpha, ncol(data)), item_params[-which(item_params == alpha)])
    names(item_params_estep) <- c(paste0("alpha", 1:ncol(data)), 
                                  names(item_params[-which(item_params == alpha)]))
  } else if (!is.null(fix_alphas)) {
    item_params_estep <- c(fix_alphas, item_params)
    names(item_params_estep) <- c(paste0("alpha", 1:ncol(data)), names(item_params))
  } else {
    item_params_estep <- item_params
    names(item_params_estep) <- names(item_params)
  }
  
  # compute post probs of the thetas with the final version of the item parameters
  post_probs <- estep_poisson_with_cov(
    data = data,
    item_params = item_params_estep,
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on
  )
  # output: a matrix with N rows (= no. of persons) and K columns (= no. of nodes)
  theta_hat <- apply(post_probs, 1, function(y) {sum(y*weights_and_nodes$x)})
  se_theta_hat <- outer(weights_and_nodes$x, theta_hat, "-")
  se_theta_hat <- t(se_theta_hat)^2 * post_probs
  se_theta_hat <- sqrt(apply(se_theta_hat, 1, sum))
  index_lower <- apply(post_probs, 1, function(x) {
    ind <- which(cumsum(x) <= .025)
    out <- ifelse(length(ind) == 0, 1, max(ind))
    return(out)
  })
  CI_lower <- weights_and_nodes$x[index_lower]
  index_upper <- apply(post_probs, 1, function(x) {min(which(cumsum(x) >= .975))})
  CI_upper <- weights_and_nodes$x[index_upper]  
  out <- data.frame(
    theta_hat = theta_hat,
    se_theta_hat = se_theta_hat,
    CI_lower = CI_lower,
    CI_upper = CI_upper
  )
  return(out)
}