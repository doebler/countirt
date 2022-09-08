
# e_step_poisson --------------------------------------------------------------------

e_step_poisson <- function(data, item_params, weights_and_nodes, 
                           item_offset = NULL, person_offset = NULL) {
  data <- as.matrix(data)
  n_items <- ncol(data)
  n_persons <- nrow(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, n_items)
  }
  
  if (is.null(person_offset)) {
    person_offset <- rep(0, n_persons)
  }
  
  PPs <- matrix(
    log(weights_and_nodes$w),
    nrow = nrow(data),
    ncol = length(weights_and_nodes$x),
    byrow = TRUE
    )
  
  for (j in 1:ncol(data)) {
    # lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]) <- old version
    log_lambdas <- alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]
    log_lambdas <- outer(person_offset, log_lambdas, "+") # result is a N * K matrix
    # if we don't have person offsets just with a lot of repetition
    lambdas <- exp(log_lambdas)
    # PPs <- PPs + outer(data[,j], lambdas, dpois, log = TRUE) <- old version
    PPs <- PPs + apply(lambdas, 2, function(x){dpois(data[,j], x, log=TRUE)})
  }
  
  PPs <- exp(PPs)
  
  PPs <- PPs / rowSums(PPs)
  
  # output should be a matrix with N rows and K cols
  return(PPs)
}

# grad_poisson -----------------------------------------------------------------------

grad_poisson <- function(item_params, PPs, weights_and_nodes, data, 
                         item_offset = NULL, person_offset = NULL) {
  data <- as.matrix(data)
  n_items <- ncol(data)
  n_persons <- nrow(data)
  alphas <- item_params[grepl("alpha", names(item_params))]
  deltas <- item_params[grepl("delta", names(item_params))]
  grad_alphas <- numeric(length(alphas))
  grad_deltas <- numeric(length(deltas))
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, n_items)
  }
  
  if (is.null(person_offset)) {
    person_offset <- rep(0, n_persons)
  }
  
  for (j in 1:ncol(data)) {
    # lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]) <- old version
    log_lambdas <- alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]
    log_lambdas <- outer(person_offset, log_lambdas, "+") # result is a N * K matrix
    # if we don't have person offsets just with a lot of repetition
    lambdas <- exp(log_lambdas)
    #x_minus_lambda <- outer(data[,j], lambdas, "-") <- old version
    x_minus_lambda <- apply(lambdas, 2, function(x){data[,j]-x})
    matrix_nodes <- matrix(
      weights_and_nodes$x,
      nrow = nrow(data), 
      ncol = length(weights_and_nodes$x),
      byrow = TRUE
      )
    x_minus_lambda_times_pp <- x_minus_lambda * PPs
    grad_alphas[j] <- sum(matrix_nodes*x_minus_lambda_times_pp)
    grad_deltas[j] <- sum(x_minus_lambda_times_pp)
  }
  
  out <- c(grad_alphas, grad_deltas)
  return(out)
}

# grad_poisson_fixalphas --------------------------------------------------------------

grad_poisson_fixalphas <- function(item_params, PPs, weights_and_nodes, 
                                   data, fix_alphas, 
                                   item_offset = NULL, person_offset = NULL) {
  data <- as.matrix(data)
  n_items <- ncol(data)
  n_persons <- nrow(data)
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params))]
  grad_deltas <- numeric(length(deltas))
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, n_items)
  }
  
  if (is.null(person_offset)) {
    person_offset <- rep(0, n_persons)
  }
  
  for (j in 1:ncol(data)) {
    # lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]) <- old version
    log_lambdas <- alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]
    log_lambdas <- outer(person_offset, log_lambdas, "+") # result is a N * K matrix
    # if we don't have person offsets just with a lot of repetition
    lambdas <- exp(log_lambdas)
    #x_minus_lambda <- outer(data[,j], lambdas, "-") <- old version
    x_minus_lambda <- apply(lambdas, 2, function(x){data[,j]-x})
    matrix_nodes <- matrix(
      weights_and_nodes$x,
      nrow = nrow(data), 
      ncol = length(weights_and_nodes$x),
      byrow = TRUE
    )
    x_minus_lambda_times_pp <- x_minus_lambda * PPs
    grad_deltas[j] <- sum(x_minus_lambda_times_pp)
  }
  
  return(grad_deltas)
}

# grad_poisson_samealpha --------------------------------------------------------------

grad_poisson_samealpha <- function(item_params, PPs, weights_and_nodes, 
                                   data, 
                                   item_offset = NULL, person_offset = NULL) {
  data <- as.matrix(data)
  n_items <- ncol(data)
  n_persons <- nrow(data)
  alphas <- rep(item_params[grepl("alpha", names(item_params))], ncol(data))
  deltas <- item_params[grepl("delta", names(item_params))]
  grad_deltas <- numeric(length(deltas))
  grad_alpha <- 0
  
  if (is.null(item_offset)) {
    item_offset <- rep(0, n_items)
  }
  
  if (is.null(person_offset)) {
    person_offset <- rep(0, n_persons)
  }
  
  for (j in 1:ncol(data)) {
    # lambdas <- exp(alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]) <- old version
    log_lambdas <- alphas[j] * weights_and_nodes$x + deltas[j] + item_offset[j]
    log_lambdas <- outer(person_offset, log_lambdas, "+") # result is a N * K matrix
    # if we don't have person offsets just with a lot of repetition
    lambdas <- exp(log_lambdas)
    #x_minus_lambda <- outer(data[,j], lambdas, "-") <- old version
    x_minus_lambda <- apply(lambdas, 2, function(x){data[,j]-x})
    matrix_nodes <- matrix(
      weights_and_nodes$x,
      nrow = nrow(data), 
      ncol = length(weights_and_nodes$x),
      byrow = TRUE
    )
    x_minus_lambda_times_pp <- x_minus_lambda * PPs
    grad_deltas[j] <- sum(x_minus_lambda_times_pp)
    grad_alpha <- grad_alpha + sum(matrix_nodes*x_minus_lambda_times_pp)
  }
  
  out <- c(grad_alpha, grad_deltas)
  return(out)
}

# em_cycle_poisson -------------------------------------------------------------------

em_cycle_poisson <- function(data, item_params, weights_and_nodes,
                     fix_alphas = NULL, same_alpha = FALSE,
                     item_offset = NULL, person_offset = NULL,
                     ctol_maxstep = 1e-8) {
    if (!is.null(fix_alphas)) {
      # fix alphas to the provided values
      # e step
      item_params_fixa <- c(fix_alphas, item_params)
      names(item_params_fixa) <- c(paste0("alpha", 1:ncol(data)), names(item_params))
      PPs <- e_step_poisson(data, item_params_fixa, weights_and_nodes, 
                            item_offset = item_offset,
                            person_offset = person_offset)
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson_fixalphas,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        fix_alphas = fix_alphas,
        item_offset = item_offset,
        person_offset = person_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_alpha) {
      # fit the model with estimating one same alpha for all item
      # e step
      alpha <- item_params[grepl("alpha", names(item_params))]
      item_params_samea <- c(rep(alpha, ncol(data)), item_params[-alpha])
      names(item_params_samea) <- c(paste0("alpha", 1:ncol(data)), 
                                   names(item_params[-alpha]))
      PPs <- e_step_poisson(data, item_params_samea, weights_and_nodes, 
                            item_offset = item_offset,
                            person_offset = person_offset)
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson_samealpha,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset,
        person_offset = person_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    } else {
      # fit a full two parameter model
      # e step
      PPs <- e_step_poisson(data, item_params, weights_and_nodes, 
                            item_offset = item_offset,
                            person_offset = person_offset)
      
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_poisson,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        item_offset = item_offset,
        person_offset = person_offset,
        control = list(xtol = ctol_maxstep)
      )$x
    }

  return(new_item_params)
}

# run_em_poisson ----------------------------------------------------------------------


run_em_poisson <- function(init_params, data, n_nodes, thres = Inf, prob = 0,
                           maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                           convcrit = "marglik",
                           fix_alphas = NULL, same_alpha = FALSE,
                           item_offset = NULL, person_offset = NULL) {
  # item_offset should be a vector at this point of length M
  # person_offset should be a vector at this point of length N
  
  # get nodes and weights for GH quadrature
  weights_and_nodes <- quad_rule(n_nodes, thres = thres,prob = prob)
  
  new_params <- init_params
  conv <- FALSE
  iter <- 1
  
  new_ll <- 0
  marg_lls <- c()
  
  print("Start estimation...")
  
  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    old_params <- new_params
    new_params <- em_cycle_poisson(
      data, old_params, weights_and_nodes,
      ctol_maxstep = ctol_maxstep,
      fix_alphas = fix_alphas, 
      same_alpha = same_alpha,
      item_offset = item_offset,
      person_offset = person_offset
    )
    
    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <- marg_ll2(
        data = as.matrix(data), 
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        family = "poisson",
        fix_alphas = fix_alphas, 
        same_alphas = same_alpha,
        item_offset = item_offset,
        person_offset = person_offset
        )
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll2(
        data = as.matrix(data), 
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        family = "poisson",
        fix_alphas = fix_alphas, 
        same_alphas = same_alpha,
        item_offset = item_offset,
        person_offset = person_offset
        )
      marg_lls[iter] <- marg_ll
      #plot(marg_lls)
      #print(marg_lls)
    }
    
    iter <- iter + 1
  }
  
  print("Done!")
  
  out <- list(
    params = new_params,
    item_offset = item_offset,
    person_offset = person_offset,
    constraints = list(fix_alphas = fix_alphas, same_alpha = same_alpha),
    iter = iter,
    conv = conv,
    marg_ll = marg_lls
  )
  return(out)
  
}

# get_start_values_pois -----------------------------------------------------------------

get_start_values_pois <- function(data, same_alpha = FALSE, fix_alphas = NULL, 
                                  item_offset = NULL,
                                  person_offset = NULL) {
  # item_offset should be a vector at this point of length M
  # person_offset should be a vector at this point of length N
  
  init_deltas <- log(apply(data, 2, mean))
  
  if (!is.null(item_offset)) {
    init_deltas <- init_deltas - item_offset
  }
  
  if (!is.null(person_offset)) {
    init_deltas <- init_deltas - mean(person_offset)
  }
  
  if (same_alpha) {
    # just one alpha for all items
    init_alphas <- c()
    for (i in 1:ncol(data)) {
      init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
    }
    init_alphas <- mean(init_alphas)
    
    start_values <- c(init_alphas, init_deltas)
    names(start_values) <- c(
      paste0("alpha", 1:length(init_alphas)),
      paste0("delta", 1:length(init_deltas))
    )
  } else if (!is.null(fix_alphas)) { 
    # we fix alphas; for current implementation of run_em_poisson, that means
    # we just don't need them to be in the item parameter vector then, because
    # they will be supplied via the fix_alphas arguments
    # so just return the delta values here
    start_values <- init_deltas
    names(start_values) <- paste0("delta", 1:length(init_deltas))
  } else {
    # different alpha for each item (full 2PPCM)
    init_alphas <- c()
    for (i in 1:ncol(data)) {
      init_alphas[i] <- cor(data[,i], apply(data[,-i], 1, mean))
    }
    
    start_values <- c(init_alphas, init_deltas)
    names(start_values) <- c(
      paste0("alpha", 1:length(init_alphas)),
      paste0("delta", 1:length(init_deltas))
    )
  }
  
  return(start_values)
}







