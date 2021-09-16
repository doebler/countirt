library(Rcpp)
library(fastGHQuad)
library(RcppGSL)
library(RcppEigen)
library(rlist)
library(rootSolve)
library(nleqslv)
source("newem.R")
source("standard_errors.R")
sourceCpp("cmp.cpp")
sourceCpp("cmp_sim_functions.cpp")
#source("BayesComp_dcomp.R")
#source("BayesComp_lambdaZ.R")
#source("BayesComp_var.R")
#load("BayesComp_logLambda.RData")
# load("BayesComp_logZ.RData")
# load("BayesComp_var.RData")
#load("grid_lambda.rda")
load("grid_log_lambda.rda")
load("grid_log_Z.rda")
load("grid_cmp_var.rda")
load("grid_W.rda")
load("grid_R.rda")

grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
grid_log_lambda_long <- as.vector(grid_log_lambda)
grid_logZ_long <- as.vector(grid_log_Z)
grid_cmp_var_long <- as.vector(grid_cmp_var)
grid_W_long <- as.vector(grid_W)
grid_R_long <- as.vector(grid_R)
grid_W_long_nona <- ifelse(is.na(grid_W_long), 
                           max(grid_W_long, na.rm = TRUE), 
                           grid_W_long)

grid_R_long_nona <- ifelse(is.na(grid_R_long), 
                           max(grid_R_long, na.rm = TRUE), 
                           grid_R_long)

rcmp <- Vectorize(rCOM_poisson)

log_lambda_from_grid <- function(mu, nu){
  
  # for numerical stability, set mu > 500 to 500
  # Warning: only works if in comparison to mean, observation is negligably small
  # and density will be 0 anyways; dcmp is not a stable implementation of density
  # if mean > 500 and density is actually relevant and different from 0
  mu <- ifelse(mu > 200, 200, mu)
  
  if (length(nu) == 1) {
    nu_filled <- rep(nu, length(mu))
  } else {
    nu_filled <- nu
  }
  
  log_lambda <- interp_from_grid_v(grid_mus, grid_nus, grid_log_lambda_long,
                                   mu, nu_filled)
  
  return(log_lambda)
}

lambda_from_grid <- function(mu, nu){
  
  # for numerical stability, set mu > 500 to 500
  # Warning: only works if in comparison to mean, observation is negligably small
  # and density will be 0 anyways; dcmp is not a stable implementation of density
  # if mean > 500 and density is actually relevant and different from 0
  mu <- ifelse(mu > 200, 200, mu)
  
  lambda_grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  lambda_grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  
  grid_log_lambda_long <- as.vector(grid_log_lambda)
  
  if (length(nu) == 1) {
    nu_filled <- rep(nu, length(mu))
  } else {
    nu_filled <- nu
  }
  
  log_lambda <- interp_from_grid_v(lambda_grid_mus, lambda_grid_nus, grid_log_lambda_long,
                                   mu, nu_filled)
  
  return(exp(log_lambda))
}

logZ_from_grid <- function(mu, nu){
  
  mu <- ifelse(mu > 200, 200, mu)
  
  grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  
  grid_logZ_long <- as.vector(grid_log_Z)
  
  if (length(nu) == 1) {
    nu_filled <- rep(nu, length(mu))
  } else {
    nu_filled <- nu
  }
  
  logZ <- interp_from_grid_v(grid_mus, grid_nus, grid_logZ_long, mu, nu_filled)
  
  return(logZ)
}

get_var_cmp <- function(mu, nu){
  # outputs bicubic interpolant of pre-computed variance
  # mu and nu can be vectors of the same length, or
  # if nu is a scalar, then replicate it to the same length as mu
  if(length(nu)==1) nu = rep(nu, length(mu))
  
  # current fineness of the gridded pre-computed values
  mu <- ifelse(mu > 200, 200, mu)
  
  grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  
  grid_cmp_var_long <- as.vector(grid_cmp_var)
  
  # perform bicubic interpolant on logLambda and logZ
  var <- interp_from_grid_v(grid_mus, grid_nus, grid_cmp_var_long, mu, nu)
  return(var)
}


dcmp <- function(data, mu, nu, logprob = FALSE) {
  # only works for length(data) = length(mu) = length(nu)
  
  # for numerical stability, set mu > 500 to 500
  # Warning: only works if in comparison to mean, observation is negligably small
  # and density will be 0 anyways; dcmp is not a stable implementation of density
  # if mean > 500 and density is actually relevant and different from 0
  mu <- ifelse(mu > 200, 200, mu)
  
  grid_mus <- c(1e-100, seq(0.001, 1, by = 0.001), as.numeric(2:200))
  grid_nus <- c(1e-100, seq(0.01, 1, by = 0.01), seq(1.1,50,0.1))
  grid_log_lambda_long <- as.vector(grid_log_lambda)
  grid_logZ_long <- as.vector(grid_log_Z)
  
  out <- dcmp_cpp(data, mu, nu, logprob,
                  grid_mus, grid_nus, grid_mus, grid_nus,
                  grid_log_lambda_long, grid_logZ_long)
  
  return(out)
}


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

# @param resp_to_item: response vector to one item j (of all N participants), summing over
#                      this vector = summing over all participants
# @param post_prob_k: posterior probabilities for one node k for all participants (different
#                     in the numerator for each participants so we have a vector over all 
#                     participants)
# @return: expected r_jk
comp_r_jk <- function(resp_to_item_j, post_prob_k) {
  r_jk <- sum(resp_to_item_j * post_prob_k)
  return(r_jk)
}

# @param post_prob_k: posterior probabilities for one node k for all participants (different
#                     in the numerator for each participants so we have a vector over all 
#                     participants)
# @return: expected f_jk (for one node)
comp_f_jk <- function(post_prob_k) {
  f_jk <- sum(post_prob_k)
  return(f_jk)
}

# @param resp_to_item: response vector to one item j (of all N participants), summing over
#                      this vector = summing over all participants
# @param post_prob_k: posterior probabilities for one node k for all participants (different
#                     in the numerator for each participants so we have a vector over all 
#                     participants)
# @return: expected h_jk
comp_h_jk <- function(resp_to_item_j, post_prob_k) {
  h_jk <- sum(vlogFactorial(resp_to_item_j) * post_prob_k)
  return(h_jk)
}


# a vectorized version of computeW
vcomputeW <- function(lambda, mu, nu) {
  out <- numeric(length(lambda))
  if (length(mu) == length(nu)) {
    for (i in 1:length(lambda)) {
      out[i] <- computeW(lambda[i], mu[i], nu[i])
    }
  } else if (length(nu) == 1) {
    for (i in 1:length(lambda)) {
      out[i] <- out[i] <- computeW(lambda[i], mu[i], nu)
    }
  } else {
    stop("vcomputeW error: All vectors entered must be of the same length, or lambda and nu must be of the same length while length of nu is 1.")
  }
  return(out)
}
vcomputeW <- Vectorize(vcomputeW)


run_e_step <- function(data, item_params, weights_and_nodes, family, fix_disps = NULL) {
  n_items <- ncol(data)
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
  e_values <- vector(mode = "list", length = n_items)
  for(j in 1:n_items) {
    lambda_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j]) 
    # yields as many lambdas as we have nodes
    
    # compute marg probability for each person (wherer for each person, we sum over all nodes)
    # TODO put this into a function that we vectorize so we don't have to loop over persons
    marg_prob <- numeric(nrow(data))
    K <- length(lambda_j)
    N <- nrow(data)
    lambda_j_matrix <- matrix(rep(lambda_j, each = N), 
                              nrow = N, ncol = K, byrow = FALSE)
    # now when we use apply columnwise we have as many lambdas from the matrix as
    # we have responses to item j, so we can feed the response vector into the dcmp
    # function and only need to rep dispersions for length of response vector to item j
    if (family == "poisson") {
      dens_i_eachnode <- apply(lambda_j_matrix, 2, function(x){dpois(data[ , j], x)})
    } else if (family == "cmp") {
      dens_i_eachnode <- apply(lambda_j_matrix, 2, function(x){dcmp(
        as.numeric(data[ , j]), 
        x,
        rep(disps[j], N)
      )})
    }
    # in dens_i_eachnode columns represent one node each; (rows = no. of responses)
    # so here we want to get a person- (but not node-)specific marginal likelihood
    # we therefore need to go rowwise so we still have one value for each person
    # and need to sum across all nodes = all columns
    marg_prob <- apply(dens_i_eachnode, 1, function(x){sum(x*weights_and_nodes$w)})
    
    # for(i in 1:nrow(data)) {
    #   if (family == "poisson") {
    #     dens_i_eachnode <- dpois(data[i,j], lambda_j)
    #   } else if (family == "cmp") {
    #     dens_i_eachnode <- dcmp(
    #       as.numeric(rep(data[i,j], length(lambda_j))),
    #       as.numeric(lambda_j),
    #       as.numeric(rep(disps[j], length(lambda_j)))
    #     )
    #   }
    #   # so we get the density under each different nodes (as contained in lambda_j)
    #   # for the one person i and the one item j we are currently looking at
    #   marg_prob[i] <- sum(dens_i_eachnode * weights_and_nodes$w)
    #   # we then sum over nodes for one specific person (because the marginal likelihood
    #   # is person specific, but not node-specific)
    # }
    
    # lambda_j_matrix <- matrix(rep(lambda_j, each = N), 
    #                           nrow = N, ncol = K, byrow = FALSE)
    # # FIXME dimenisons can't be right here, data[,j] is not a scalar
    # if (family == "poisson") {
    #   dens <- apply(lambda_j_matrix, 2, function(x){dpois(data[ , j], x)}) 
    # } else if (family == "cmp") {
    #   dens <- apply(lambda_j_matrix, 2, function(x){dcmp(
    #     as.numeric(rep(data[ , j], length(lambda_j))), 
    #     x,
    #     rep(disps[j], length(lambda_j))
    #     )})
    # }
    
    # compute posterior probabilties for all nodes and persons (for item j) in
    # the shape of a matrix with dimensions persons (rows) + nodes (columns)
    post_probs <- computepp_allnodes(dens_i_eachnode, weights_and_nodes$w, marg_prob)
    r_j <- computer_allnodes(data[ , j], post_probs)
    f_j <- computef_allnodes(post_probs)
    if (family == "cmp") {
      h_j <- computeh_allnodes(data[ , j], post_probs)
    } else {
      h_j <- NULL
    }
    
    e_values[[j]] <- list(post_probs = post_probs, r_j = r_j, f_j = f_j, h_j = h_j)
    names(e_values)[j] <- paste0("item", j)
    # each entry in the list is for one item and that has three vectors with as many
    # elements as there are nodes
  }
  return(e_values)
}

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

# grad_e_ll_cmp <- function(item_params, e_values, weights_and_nodes) {
#   alphas <- item_params[grepl("alpha", names(item_params))]
#   deltas <- item_params[grepl("delta", names(item_params))]
#   disps <- item_params[grepl("disp", names(item_params))]
#   n_items <- length(alphas)
#   
#   grads <- numeric(length(item_params))
#   for(j in 1:length(item_params)) {
#     # for all nodes compute cmp mean
#     if(grepl("alpha", names(item_params)[j])) {
#       mu_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
#       V_j <- BayesComp_var(mu_j, disps[j])
#       # compute lambdas for current item parameter and all nodes
#       # we need the expected counts for the current item j at each of node k
#       # lambda_j is a vector with k elements here, as is r_j and f_j
#       grads[j] <- sum((weights_and_nodes$x * mu_j / V_j) * (e_values[[j]][["r_j"]] -
#                                         - mu_j * e_values[[j]][["f_j"]]))
#     } else if (grepl("delta", names(item_params)[j])){
#       mu_j <- exp(alphas[j - n_items] * weights_and_nodes$x + deltas[j - n_items])
#       V_j <- BayesComp_var(mu_j, disps[j - n_items])
#       # compute lambdas for current item parameter and all nodes
#       grads[j] <- sum((mu_j / V_j) * (e_values[[j - n_items]][["r_j"]] -
#                                          - mu_j * e_values[[j - n_items]][["f_j"]]))
#     } else {
#       mu_j <- exp(alphas[j - 2*n_items] * weights_and_nodes$x + deltas[j - 2*n_items])
#       V_j <- BayesComp_var(mu_j, disps[j - 2*n_items])
#       # dispersion parameters
#       # FIXME new version of BayesComp_lambdaZ returns lambda not logLambda
#       lambda <- exp(BayesComp_lambdaZ(mu_j, disps[j - 2*n_items])$logLambda)
#       e_logfac <- tmbElogFactorial(lambda, disps[j - 2*n_items])
#       # to avoid getting NAs at the borders which will cause gradient to fail,
#       # we will set them to 0 when they will be multiplied with a 0 in the gradient anyways
#       # (or with a values as good as zero, i.e., 2e-52, double precision)
#       e_logfac <- ifelse(
#         e_values[[j - 2*n_items]][["f_j"]] <= 2e-52,
#         ifelse(
#           e_values[[j - 2*n_items]][["f_j"]] == 0,
#           0, 
#           2e-52
#         ),
#         e_logfac
#       )
#       e_logfac_x <- tmbYlogFactorial(lambda, disps[j - 2*n_items])
#       e_logfac_x <- ifelse(
#         (e_values[[j - 2*n_items]][["r_j"]] -
#           mu_j * e_values[[j - 2*n_items]][["f_j"]]) <= 2e-52,
#         ifelse(
#           (e_values[[j - 2*n_items]][["r_j"]] -
#              mu_j * e_values[[j - 2*n_items]][["f_j"]]) == 0,
#           0,
#           2e-52
#         ),
#         e_logfac_x
#       )
#       grads[j] <- sum((e_logfac * e_values[[j - 2*n_items]][["f_j"]]) - 
#                         e_values[[j - 2*n_items]][["h_j"]] +
#                         ((e_logfac_x - mu_j * e_logfac) / V_j) * 
#                         (e_values[[j - 2*n_items]][["r_j"]] -
#                         mu_j * e_values[[j - 2*n_items]][["f_j"]]))
#     }
#   }
#   print(grads)
#   return(grads)
# }

# we estimate the log dispersion instead of the dispersion because we need the
# to bound the parameter space 
grad_e_ll_cmp_logdisp <- function(item_params, e_values, weights_and_nodes) {
  # print(item_params)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- length(alphas)
  n_nodes <- length(weights_and_nodes$x)
  
  grads <- numeric(length(item_params))
  for(j in 1:length(item_params)) {
    # for all nodes compute cmp mean
    if(grepl("alpha", names(item_params)[j])) {
      mu_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
      #V_j <- BayesComp_var(mu_j, disps[j])
      V_j <- get_var_cmp(mu_j, disps[j])
      # compute lambdas for current item parameter and all nodes
      # we need the expected counts for the current item j at each of node k
      # lambda_j is a vector with k elements here, as is r_j and f_j
      grads[j] <- sum((weights_and_nodes$x * mu_j / V_j) * (e_values[[j]][["r_j"]] -
                                                              mu_j * e_values[[j]][["f_j"]]))
    } else if (grepl("delta", names(item_params)[j])){
      mu_j <- exp(alphas[j - n_items] * weights_and_nodes$x + deltas[j - n_items])
      V_j <- get_var_cmp(mu_j, disps[j - n_items])
      # compute lambdas for current item parameter and all nodes
      grads[j] <- sum((mu_j / V_j) * (e_values[[j - n_items]][["r_j"]] -
                                        mu_j * e_values[[j - n_items]][["f_j"]]))
    } else {
      # dispersion parameters
      mu_j <- exp(alphas[j - 2*n_items] * weights_and_nodes$x + deltas[j - 2*n_items])
      V_j <- get_var_cmp(mu_j, disps[j - 2*n_items])
      lambda <- lambda_from_grid(mu_j, disps[j - 2*n_items])
      W_j <- rep(NA, length(n_nodes))
      indexW <- (e_values[[j - 2*n_items]][["r_j"]] > 1e-8) | 
        ((mu_j * e_values[[j - 2*n_items]][["f_j"]]) > 1e-8)
      W_j[indexW] <- vcomputeW(lambda[indexW], mu_j[indexW], disps[j - 2*n_items])
      # whenever W is NA, this should coincide with posterior porbabilities of next to 0
      # so that we check how large r and f are and whenever they are 0 or next to 0, we set 
      # to either a very small value or 0
      frac_r_W <- rep(NA, length(n_nodes))
      index1 <- (e_values[[j - 2*n_items]][["r_j"]] <= 1e-8)
      index2 <- (e_values[[j - 2*n_items]][["r_j"]] == 0)
      index3 <- (e_values[[j - 2*n_items]][["r_j"]] > 1e-8)
      frac_r_W[index1] <- 1e-8
      frac_r_W[index2] <- 0
      frac_r_W[index3] <- e_values[[j - 2*n_items]][["r_j"]][index3] / W_j[index3]
      # frac_r_W <- ifelse(
      #   e_values[[j - 2*n_items]][["r_j"]] <= 1e-30,
      #   ifelse(
      #     e_values[[j - 2*n_items]][["r_j"]] == 0,
      #     0, 
      #     1e-30
      #   ),
      #   e_values[[j - 2*n_items]][["r_j"]] / W_j
      # )
      frac_muf_W <- rep(NA, length(n_nodes))
      index3 <- ((mu_j * e_values[[j - 2*n_items]][["f_j"]]) <= 1e-8)
      index4 <- ((mu_j * e_values[[j - 2*n_items]][["f_j"]]) == 0)
      index5 <- ((mu_j * e_values[[j - 2*n_items]][["f_j"]]) > 1e-8)
      frac_muf_W[index3] <- 1e-8
      frac_muf_W[index4] <- 0
      frac_muf_W[index5] <- (mu_j[index5] * e_values[[j - 2*n_items]][["f_j"]][index5]) / 
        W_j[index5]
      # frac_muf_W <- ifelse(
      #   (mu_j * e_values[[j - 2*n_items]][["f_j"]]) <= 1e-30,
      #   ifelse(
      #     (mu_j * e_values[[j - 2*n_items]][["f_j"]]) == 0,
      #     0, 
      #     1e-30
      #   ),
      #   (mu_j * e_values[[j - 2*n_items]][["f_j"]]) / W_j
      # )
      
      # same logic applies as for W's, just make sure we set cases to 0 where logfac function fails
      # but posterior probabilities are as good as 0 anyways
      e_logfac <- rep(NA, length(n_nodes))
      index6 <- (e_values[[j - 2*n_items]][["f_j"]] <= 1e-8)
      index7 <- (e_values[[j - 2*n_items]][["f_j"]] == 0)
      index8 <- (e_values[[j - 2*n_items]][["f_j"]] > 1e-8)
      e_logfac[index6] <- 1e-8
      e_logfac[index7] <- 0
      e_logfac[index8] <- tmbElogFactorial(lambda[index8], disps[j - 2*n_items])
      # e_logfac <- ifelse(
      #   e_values[[j - 2*n_items]][["f_j"]] <= 1e-30,
      #   ifelse(
      #     e_values[[j - 2*n_items]][["f_j"]] == 0,
      #     0, 
      #     1e-30
      #   ),
      #   tmbElogFactorial(lambda, disps[j - 2*n_items])
      # )
      
      grads[j] <- sum(disps[j - 2*n_items] * 
                        (frac_r_W - frac_muf_W - e_values[[j - 2*n_items]][["h_j"]] + 
                           e_logfac * e_values[[j - 2*n_items]][["f_j"]]))
    }
  }
  #  print(grads)
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

grad_e_ll_cmp_fixdisps <- function(item_params, e_values, weights_and_nodes, fix_disps) {
  #print("item params:")
  #print(item_params)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  disps <- fix_disps
  n_items <- length(alphas)
  
  grads <- numeric(length(item_params))
  for(j in 1:length(item_params)) {
    # for all nodes compute cmp mean
    if(grepl("alpha", names(item_params)[j])) {
      mu_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
      V_j <- get_var_cmp(mu_j, disps[j])
      # compute lambdas for current item parameter and all nodes
      # we need the expected counts for the current item j at each of node k
      # lambda_j is a vector with k elements here, as is r_j and f_j
      grads[j] <- sum((weights_and_nodes$x * mu_j / V_j) * (e_values[[j]][["r_j"]] -
                                                              mu_j * e_values[[j]][["f_j"]]))
    } else if (grepl("delta", names(item_params)[j])){
      mu_j <- exp(alphas[j - n_items] * weights_and_nodes$x + deltas[j - n_items])
      V_j <- get_var_cmp(mu_j, disps[j - n_items])
      # compute lambdas for current item parameter and all nodes
      grads[j] <- sum((mu_j / V_j) * (e_values[[j - n_items]][["r_j"]] -
                                        mu_j * e_values[[j - n_items]][["f_j"]]))
    } 
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

expect_ll_cmp <- function(item_params, e_values, weights_and_nodes) {
  # print(item_params)
  # s <<- s + 1
  # print(s)
  #e_values[[j - 2*n_items]][["h_j"]]
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  log_disps <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(log_disps)
  n_items <- length(alphas)
  out <- 0
  
  for (j in 1:n_items) {
    mu_j <- exp(alphas[j] * weights_and_nodes$x + deltas[j])
    log_lambda_j <- log(lambda_from_grid(mu_j, disps[j]))
    log_Z <- logZ_from_grid(mu_j, disps[j])
    out <- out + sum(
      e_values[[j]][["r_j"]] * log_lambda_j -
        disps[j] * e_values[[j]][["h_j"]] -
        e_values[[j]][["f_j"]] * log_Z
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

# for convergence check evaluate observed likelihood instead of expected likelihood
# because it is easier to evaluate (even though harder to maximize)
# and it is optimized by the same parameters
marg_ll <- function(data, item_params, weights_and_nodes, family, fix_disps = NULL) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  if (family == "poisson") {
    if (length(item_params) == 2*ncol(data)) {
      alphas <- item_params[grepl("alpha", names(item_params))]
    } else {
      alphas <- rep(1, ncol(data))
    }
  } else if (family == "cmp") {
    alphas <- item_params[grepl("alpha", names(item_params))]
    if (is.null(fix_disps)) {
      log_disps <- item_params[grepl("log_disp", names(item_params))]
      disps <- exp(log_disps)
    } else {
      disps <- fix_disps
    }
  }
  
  
  if (family == "poisson") {
    # function to compute integral with quadrature over
    f <- function(z, data, alphas, deltas) {
      out <- 0
      for (j in 1:n_items) {
        lambda <- exp(alphas[j] * z + deltas[j])
        out <- out + (dpois(data[,j], lambda, log = TRUE))
      }
      return(out)
    }
    
    log_marg_prob <- numeric(n_persons)
    for (i in 1:n_persons) {
      log_marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                                 data = data[i, , drop = FALSE], 
                                 alphas = alphas, deltas = deltas)
    }
  } else if (family == "cmp") {
    # function to compute integral with quadrature over
    f <- function(z, data, alphas, deltas, disps) {
      out <- 0
      for (j in 1:n_items) {
        mu <- exp(alphas[j] * z + deltas[j])
        out <- out + (dcmp(rep(as.numeric(data[,j]), length(mu)),
                           mu, rep(disps[j], length(mu)), 
                           logprob = TRUE))
        #out <- out * (BayesComp_dcomp(data[,j], mu, disps[j]))
      }
      return(out)
    }
    
    log_marg_prob <- numeric(n_persons)
    for (i in 1:n_persons) {
      log_marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                                 data = data[i, , drop = FALSE], 
                                 alphas = alphas, deltas = deltas, 
                                 disps = disps)
    }
  }
  
  ll <- sum(log_marg_prob)
  return(ll)
}

marg_ll2 <- function(data, item_params, weights_and_nodes, family, fix_disps = NULL,
                     fix_alphas = NULL, interp_method = "bicubic",
                     same_disps = FALSE, same_alphas = FALSE) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params))]
  if (family == "poisson") {
    if (is.null(fix_alphas)) {
      if (same_alphas) {
        alpha <- item_params[grepl("alpha", names(item_params))]
        alphas <- rep(alpha, n_items) 
      } else {
        alphas <- item_params[grepl("alpha", names(item_params))]
      }
    } else {
      alphas <- fix_alphas
    }
  } else if (family == "cmp") {
    if (is.null(fix_alphas)) {
      if (same_alphas) {
        alpha <- item_params[grepl("alpha", names(item_params))]
        alphas <- rep(alpha, n_items)
      } else {
        alphas <- item_params[grepl("alpha", names(item_params))]
      }
    } else {
      alphas <- fix_alphas
    }
    if (is.null(fix_disps)) {
      if (same_disps) {
        log_disp <- item_params[grepl("log_disp", names(item_params))]
        disps <- rep(exp(log_disp), n_items)
      } else {
        log_disps <- item_params[grepl("log_disp", names(item_params))]
        disps <- exp(log_disps)
      }
    } else {
      disps <- fix_disps
    }
  }
  
  
  if (family == "poisson") {
    # function to compute integral with quadrature over
    f <- function(z, data, alphas, deltas) {
      out <- 0
      for (j in 1:n_items) {
        lambda <- exp(alphas[j] * z + deltas[j])
        out <- out + (dpois(data[,j], lambda, log = TRUE))
      }
      return(exp(out))
    }
    
    marg_prob <- numeric(n_persons)
    for (i in 1:n_persons) {
      marg_prob[i] <- ghQuad(f, rule = weights_and_nodes,
                             data = data[i, , drop = FALSE], 
                             alphas = alphas, deltas = deltas)
    }
    ll <- sum(log(marg_prob))
  } else if (family == "cmp") {
    if (interp_method == "bicubic") {
      ll <- marg_ll_cpp(data = as.matrix(data),
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
      # then interpolation method is linear
      # here we don't cap mu, so we extrapolate beyond grid values
      ll <- marg_ll_cpp_lininterp(data = as.matrix(data),
                                  alphas = alphas,
                                  deltas = deltas, 
                                  disps = disps,
                                  nodes = weights_and_nodes$x,
                                  weights = weights_and_nodes$w,
                                  grid_mus = grid_mus,  
                                  grid_nus = grid_nus, 
                                  grid_logZ_long = grid_logZ_long,
                                  grid_log_lambda_long = grid_log_lambda_long)
    }
    
  }
  return(ll)
}

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
                                  thres = Inf, prob = 0) {
  # get weights and nodes
  weights_and_nodes <- quad_rule(n_nodes, thres = thres,prob = prob)
  # compute post probs of the thetas with the final version of the item parameters
  post_probs <- newem_estep2(
    data = data,
    item_params = item_params,
    weights_and_nodes = weights_and_nodes
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

# get_start_values -----------------------------------------------------------------

get_start_values <- function(data, init_disp_one = TRUE, same_alpha = FALSE) {
  init_deltas <- log(apply(data, 2, mean))
  
  if (same_alpha) {
    # just one alpha for all items
    init_alphas <- c()
    for (i in 1:ncol(data)) {
      init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
    }
    init_alphas <- mean(init_alphas)
  } else {
    # different alpha for each item
    init_alphas <- c()
    for (i in 1:ncol(data)) {
      init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
    }
  }
  
  if (init_disp_one) {
    init_logdisps <- log(rep(1, length(init_deltas)))
  } else {
    # we could try and get start values with:
    init_logdisps <- apply(data, 2, mean) / apply(data, 2, var)
    # but actually, we get the same results upon trying out
    # with disp = 1
  }
  
  start_values <- c(init_alphas, init_deltas, init_logdisps)
  names(start_values) <- c(
    paste0("alpha", 1:length(init_alphas)),
    paste0("delta", 1:length(init_deltas)),
    paste0("log_disp", 1:length(init_logdisps))
  )
  return(start_values)
}

# gen_test_data ----------------------------------------------------------------------------------------

gen_test_data <- function(true_alphas, true_deltas, true_abilities, n_persons) {
  out <- data.frame(
    item1 = numeric(n_persons)
  )
  for (j in 2:length(true_alphas)) {
    out[[paste0("item", j)]] <- numeric(n_persons)
  }
  for (j in 1:length(true_alphas)) {
    lambdas <- exp(true_alphas[j] * true_abilities + true_deltas[j])
    out[[paste0("item", j)]] <- rpois(n_persons, lambdas)
  }
  return(out)
}

# gen_test_data_cmp ----------------------------------------------------------------------

gen_test_data_cmp <- function(true_alphas, true_deltas, true_disps, 
                              true_abilities, n_persons) {
  n_items <- length(true_alphas)
  abilities <- matrix(
    true_abilities,
    ncol = n_items,
    nrow = n_persons,
    byrow = FALSE
  )
  true_alphas_m <- matrix(
    rep(true_alphas, each = n_persons), 
    ncol = n_items,
    nrow = n_persons, 
    byrow = FALSE
  )
  true_deltas_m <- matrix(
    rep(true_deltas, each = n_persons), 
    ncol = n_items,
    nrow = n_persons, 
    byrow = FALSE
  )
  true_disps_v <- rep(true_disps, each = n_persons)
  mus <- exp(true_alphas_m * abilities + true_deltas_m)
  mus <- as.vector(mus) # long format, first all people for item 1, then 
  # all people for item 2, etc.
  lambdas <- exp(log_lambda_from_grid(mus, true_disps_v))
  # sample from CMP distribution
  out <- rcmp(runif(n_persons*n_items), lambdas, true_disps_v) 
  out <- as.data.frame(matrix(out, ncol = n_items)) 
  return(out)
}






