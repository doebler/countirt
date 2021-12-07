
# get_resp_patterns_pcov_cat -------------------------------------------------------------

# this functions returns a list with the dummy coding matrices for each categorical person
# covriate 
get_resp_patterns_pcov_cat <- function(n_levels) {
  out <- matrix(0, ncol = n_levels, nrow = n_levels)
  for (l in 2:n_levels) {
    out[l,l] <- 1
  }
  return(out)
}

# make_resp_patterns_mat ------------------------------------------------------------

make_resp_patterns_mat <- function(resp_pattern_list, n_resp_patterns, num_levels_p_cov) {
  out <- matrix(NA, ncol = sum(num_levels_p_cov), nrow = n_resp_patterns)
  n_cov <- length(resp_pattern_list)
  cov_levels <- lapply(num_levels_p_cov, function(x){1:x})
  index_combis <- expand.grid(cov_levels)
  for (l in 1:n_resp_patterns) {
    out[l,] <- do.call(
      cbind,
      lapply(1:n_cov, function(x){resp_pattern_list[[x]][index_combis[l,x],,drop=FALSE]})
    )
  }
  # remove the columns for the reference groups where all entries are 0
  out <- out[, colSums(out) > 0, drop = FALSE]
  return(out)
}







