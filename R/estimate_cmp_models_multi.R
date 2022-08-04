









# run_em_multi ---------------------------------------------------------------------------

# TODO hier weitermachen und von unten nach oben implementieren auf basis der poisson 
# funktionen
run_em_multi <- function() {
  
}

# get_start_values_multi -----------------------------------------------------------------
# TODO disp_constraints implementieren, also das wir die fixieren koennen und auch
# equality constraints einbauen koennen
get_start_values_multi <- function(data, 
                                   n_traits, 
                                   alpha_constraints = NULL,
                                   n_nodes = NULL, n_samples = NULL, 
                                   em_type = c("gh", "mc"), 
                                   fcov_prior = NULL, truncate_grid = TRUE,
                                   penalize = c("none", "ridge", "lasso"), 
                                   penalize_lambda = NULL,
                                   maxiter = 100,
                                   nsim = 1000) {
  # use maxiter here to not let the poisson model run til full convergence but just
  # enough to get some sensible start values as the poisson em will also already take
  # quite a bit of time here
  
  # for CMP start values, we fit a Poisson model and get deltas and alphas from there
  init_values_pois <- get_start_values_pois_multi(
    data = data, 
    n_traits = n_traits, 
    alpha_constraints = alpha_constraints
  )
  fit_pois <- run_em_poisson_multi(
    data = data, 
    init_params = init_values_pois, 
    n_traits = n_traits, 
    n_nodes = n_nodes, 
    n_samples = n_samples, 
    em_type = em_type, 
    fcov_prior = fcov_prior,
    truncate_grid = truncate_grid,
    penalize = penalize, 
    penalize_lambda = penalize_lambda,
    maxiter = maxiter,
    alpha_constraints = alpha_constraints
  )
  init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
  # the outputted alphas are the whole matrix, with constrained alphas fixed at 
  # their constrained values, so e.g. fixed at 0; i thus need to pass alpha_constraints
  # to all subsequent functions and check which of these values need never be
  # estimated because they were fixed at a certain values
  init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  
  # start values for log nus
  init_logdisps <- c()
  sim_abilities <- mvrnorm(nsim, rep(0, n_traits), diag(rep(1, n_traits)))
  for (i in 1:ncol(data)) {
    alphas_for_item_i <- init_alphas[grepl(paste0("alpha", i), )]
    mu <- exp(init_deltas[i] + alphas_for_item_i%*%sim_abilities)
    sim <- rpois(nsim, mu)
    init_logdisps[i] <- log((var(sim) / var(data[,i])))
  }
  
  start_values <- c(init_alphas, init_deltas, init_logdisps)
  names(start_values) <- c(
    names(init_alphas),
    names(init_deltas),
    paste0("log_disp", 1:length(init_logdisps))
  )
  return(start_values)
}
