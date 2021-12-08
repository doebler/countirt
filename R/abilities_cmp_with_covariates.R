# grad_ll_cmp_ability_1P ---------------------------------------------------------------------------------------

# first derivative of log probability of the response vector for one participant
# @param ability: ability value for the one participant we are looking at
# @param data_1P: vector of responses of one person to all items, length = no. of items
# @param item_params: item parameters estimates as returned by the model estimation; must be named
#                     alpha1, ..., alphaM, delta1, ..., deltaM, log_disp1, ..., log_dispM
# grad_ll_cmp_ability_1P <- function(ability, data_1P, item_params) {
#   n_items <- length(data_1P)
#   alphas <- item_params[grepl("alpha", names(item_params))]
#   deltas <- item_params[grepl("delta", names(item_params))]
#   log_disps <- item_params[grepl("log_disp", names(item_params))]
#   disps <- exp(log_disps)
#   
#   grad <- grad_ll_cmp_ability_1P_cpp(
#     ability = ability,
#     alphas = alphas, 
#     deltas = deltas, 
#     disps = disps, 
#     data = data_1P,
#     grid_mus = grid_mus, 
#     grid_nus = grid_nus, 
#     grid_cmp_var_long = grid_cmp_var_long,
#     grid_log_lambda_long = grid_log_lambda_long,
#     max_mu = 200,
#     min_mu = 0.001
#   )
#   
#   return(grad) 
# }

# ll_cmp_ability_1P ------------------------------------------------------------------------------------------

# log probability of the response vector for one participant
# @param ability: ability value for the one participant we are looking at
# @param data_1P: vector of responses of one person to all items, length = no. of items
# @param item_params: item parameters estimates as returned by the model estimation; must be named
#                     alpha1, ..., alphaM, delta1, ..., deltaM, log_disp1, ..., log_dispM
# ll_cmp_ability_1P <- function(ability, data_1P, item_params) {
#   n_items <- length(data_1P)
#   alphas <- item_params[grepl("alpha", names(item_params))]
#   deltas <- item_params[grepl("delta", names(item_params))]
#   log_disps <- item_params[grepl("log_disp", names(item_params))]
#   disps <- exp(log_disps)
#   
#   out <- 0
#   for (j in 1:n_items) {
#     mu <- exp(alphas[j]*ability + deltas[j])
#     out <- out + dcmp(as.numeric(data[j]), mu, disps[j], logprob = TRUE)
#   }
#   
#   return(out)
# }

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
# get_ability_params <- function(data, item_params, init_values = 0, ctol = 1e-8) {
#   n_persons <- nrow(data)
#   if (length(init_values) == 1) {
#     init_values <- rep(init_values, n_persons)
#   }
#   abilities <- numeric(n_persons)
#   se_abilities <- numeric(n_persons)
#   
#   for (i in 1:n_persons) {
#     
#     optimum <- nleqslv(
#       x = init_values[i],
#       fn = grad_ll_cmp_ability_1P,
#       data_1P = as.numeric(data[i,]), 
#       item_params = item_params,
#       control = list(xtol = ctol),
#       jacobian = TRUE
#     )
#     
#     abilities[i] <- optimum$x
#     se_abilities[i] <- sqrt(solve(-optimum$jac))
#     
#   }
#   
#   out <- data.frame(
#     abilities = abilities,
#     se_abilities = se_abilities
#   )
#   
#   return(out)
# }

# get_ability_params_pp_with_cov ------------------------------------------------------------

get_ability_params_pp_with_cov <- function(data, item_params, n_nodes,
                                           p_covariates, i_covariates,
                                           i_cov_on = c("alpha", "delta", "log_disp"),
                                           p_cov_cat = TRUE,
                                           resp_patterns_matrix = NULL,
                                           fix_disps = NULL, fix_alphas = NULL,
                                           same_disps = FALSE, same_alphas = FALSE,
                                           thres = Inf, prob = 0) {
  # get weights and nodes
  weights_and_nodes <- quad_rule(n_nodes, thres = thres, prob = prob)
  
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
      item_params_estep <- c(alphas, deltas, rep(log_disp, n_items), 
                                betas_i_alpha, betas_i_delta)
      names(item_params_estep) <- c(
        names(item_params[grepl("alpha", names(item_params))]),
        names(item_params[grepl("delta", names(item_params))]),
        paste0("log_disp", 1:n_items),
        names(item_params[grepl("beta_i_alpha", names(item_params))]),
        names(item_params[grepl("beta_i_delta", names(item_params))])
      )
    } else {
      # either we have person covariates or just item covariates on one parameter, so just beta_i
      item_params_estep <- c(alphas, deltas, rep(log_disp, n_items), betas_p, betas_i)
      names(item_params_estep) <- c(
        names(item_params[grepl("alpha", names(item_params))]),
        names(item_params[grepl("delta", names(item_params))]),
        paste0("log_disp", 1:n_items),
        names(item_params[grepl("beta_p", names(item_params))]),
        names(item_params[grepl("beta_i", names(item_params))])
      )
    }
  } else if (!is.null(fix_disps)) {
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
      item_params_estep <- c(alphas, deltas, log(fix_disps), 
                                betas_i_alpha, betas_i_delta)
      names(item_params_estep) <- c(
        names(item_params[grepl("alpha", names(item_params))]),
        names(item_params[grepl("delta", names(item_params))]),
        paste0("log_disp", 1:n_items),
        names(item_params[grepl("beta_i_alpha", names(item_params))]),
        names(item_params[grepl("beta_i_delta", names(item_params))])
      )
    } else {
      # either we have person covariates or just item covariates on one parameter, so just beta_i
      item_params_estep <- c(alphas, deltas, log(fix_disps), betas_p, betas_i)
      names(item_params_estep) <- c(
        names(item_params[grepl("alpha", names(item_params))]),
        names(item_params[grepl("delta", names(item_params))]),
        paste0("log_disp", 1:n_items),
        names(item_params[grepl("beta_p", names(item_params))]),
        names(item_params[grepl("beta_i", names(item_params))])
      )
    }
  } else {
    item_params_estep <- item_params
    names(item_params_estep) <- names(item_params)
  }
  
  # compute post probs of the thetas with the final version of the item parameters
  post_probs <- estep_cmp_with_cov(
    data = data,
    item_params = item_params,
    weights_and_nodes = weights_and_nodes,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    resp_patterns_matrix = resp_patterns_matrix
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







