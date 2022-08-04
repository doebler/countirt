














# function for setting start values
# TODO think about setting down convergence requirements (tol) on poisson
# estimation to save some time and also set down max iter
# this is especially relevent for mc, where we dont wanna run super many iters
# just due to random fluctuation

# get_start_values_multi -----------------------------------------------------------------

get_start_values_multi <- function(data, nodes = 121, nsim = 1000,
                             init_disp_one = TRUE, same_alpha = FALSE) {
  # for CMP start values, we fit a Poisson model and get deltas and alphas from there
  if (same_alpha) {
    # just one alpha for all items
    init_values_pois <- get_start_values_pois(data, same_alpha = TRUE)
    fit_pois <- run_em_poisson(data, init_values_pois, nodes, same_alpha = TRUE)
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  } else {
    # different alpha for each item
    init_values_pois <- get_start_values_pois(data)
    fit_pois <- run_em_poisson(data, init_values_pois, nodes)
    init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params))]
    init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params))]
  }
  
  init_logdisps<-c()
  sim_abilities=rnorm(nsim)
  for (i in 1:ncol(data)) {
    if (same_alpha) {
      mu <- exp(init_deltas[i] + init_alphas*sim_abilities)
    } else {
      mu <- exp(init_deltas[i] + init_alphas[i]*sim_abilities)
    }
    sim <- rpois(nsim, mu)
    init_logdisps[i] <- log((var(sim) / var(data[,i])))
  }
  
  start_values <- c(init_alphas, init_deltas, init_logdisps)
  names(start_values) <- c(
    paste0("alpha", 1:length(init_alphas)),
    paste0("delta", 1:length(init_deltas)),
    paste0("log_disp", 1:length(init_logdisps))
  )
  return(start_values)
}