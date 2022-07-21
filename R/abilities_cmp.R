# grad_ll_cmp_ability_1P ---------------------------------------------------------------------------------------

# first derivative of log probability of the response vector for one participant
# @param ability: ability value for the one participant we are looking at
# @param data_1P: vector of responses of one person to all items, length = no. of items
# @param item_params: item parameters estimates as returned by the model estimation; must be named
#                     alpha1, ..., alphaM, delta1, ..., deltaM, log_disp1, ..., log_dispM
grad_ll_cmp_ability_1P <- function(ability, data_1P, item_params) {
  n_items <- length(data_1P)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  grad <- grad_ll_cmp_ability_1P_cpp(
    ability = ability,
    alphas = alphas, 
    deltas = deltas, 
    disps = disps, 
    data = data_1P,
    grid_mus = grid_mus, 
    grid_nus = grid_nus, 
    grid_cmp_var_long = grid_cmp_var_long,
    grid_log_lambda_long = grid_log_lambda_long,
    max_mu = 200,
    min_mu = 0.001
  )
  
  return(grad) 
}

# ll_cmp_ability_1P ------------------------------------------------------------------------------------------

# log probability of the response vector for one participant
# @param ability: ability value for the one participant we are looking at
# @param data_1P: vector of responses of one person to all items, length = no. of items
# @param item_params: item parameters estimates as returned by the model estimation; must be named
#                     alpha1, ..., alphaM, delta1, ..., deltaM, log_disp1, ..., log_dispM
ll_cmp_ability_1P <- function(ability, data_1P, item_params) {
  n_items <- length(data_1P)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  
  out <- 0
  for (j in 1:n_items) {
    mu <- exp(alphas[j]*ability + deltas[j])
    out <- out + dcmp(as.numeric(data[j]), mu, disps[j], logprob = TRUE)
  }
  
  return(out)
}

# get_ability_params --------------------------------------------------------------

# @param data: data matrix with as many columns as we have items and as many rows as we have columns
# @param item_params: item parameters estimates as returned by the model estimation; must be named
#                     alpha1, ..., alphaM, delta1, ..., deltaM, log_disp1, ..., log_dispM
# @param init_values: how to initialize ability parameter estimation, defaults
#                     0, so for all persons, we start with a 0, if to be initalized
#                     with different values, please provide a vector of values
#                     for each person, if only a scalar is provided, it will be
#                     used for all persons
# @return a dataframe with two columns: one with ability parameters and one with 
#         corresponding SE's for each person
get_ability_params <- function(data, item_params, init_values = 0, ctol = 1e-8) {
  n_persons <- nrow(data)
  if (length(init_values) == 1) {
    init_values <- rep(init_values, n_persons)
  }
  abilities <- numeric(n_persons)
  se_abilities <- numeric(n_persons)
  
  for (i in 1:n_persons) {
    
    optimum <- nleqslv(
      x = init_values[i],
      fn = grad_ll_cmp_ability_1P,
      data_1P = as.numeric(data[i,]), 
      item_params = item_params,
      control = list(xtol = ctol),
      jacobian = TRUE
    )
    
    abilities[i] <- optimum$x
    se_abilities[i] <- sqrt(solve(-optimum$jac))
    
  }
  
  out <- data.frame(
    abilities = abilities,
    se_abilities = se_abilities
  )
  
  return(out)
}

# get_ability_params_pp ------------------------------------------------------------

get_ability_params_pp <- function(data, item_params, n_nodes,
                                  fix_disps = NULL, fix_alphas = NULL,
                                  same_disps = FALSE, same_alphas = FALSE,
                                  thres = Inf, prob = 0,
                                  item_offset = NULL) {
  # get weights and nodes
  weights_and_nodes <- quad_rule(n_nodes, thres = thres,prob = prob)
  
  # prep the item parameters for e step
  if (same_alphas) {
    alpha <- item_params[grepl("alpha", names(item_params))]
    item_params_estep <- c(rep(alpha, ncol(data)), item_params[-which(item_params == alpha)])
    names(item_params_estep) <- c(paste0("alpha", 1:ncol(data)), 
                                  names(item_params[-which(item_params == alpha)]))
  } else if (!is.null(fix_alphas)) {
    item_params_estep <- c(fix_alphas, item_params)
    names(item_params_estep) <- c(paste0("alpha", 1:ncol(data)), names(item_params))
  } else if (same_disps) {
    log_disp <- item_params[grepl("log_disp", names(item_params))]
    item_params_estep <- c(item_params[-which(item_params == log_disp)], 
                           rep(log_disp, ncol(data)))
    names(item_params_estep) <- c(names(item_params[-which(item_params == log_disp)]),
                                  paste0("log_disp", 1:ncol(data)))
  } else if (!is.null(fix_disps)) {
    item_params_estep <- c(item_params, log(fix_disps))
    names(item_params_estep) <- c(names(item_params), 
                                  paste0("log_disp", 1:ncol(data)))
  } else {
    item_params_estep <- item_params
    names(item_params_estep) <- names(item_params)
  }
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, ncol(data))
  }
  
  # compute post probs of the thetas with the final version of the item parameters
  post_probs <- newem_estep2(
    data = data,
    item_params = item_params_estep,
    weights_and_nodes = weights_and_nodes,
    item_offset = item_offset
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







