
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
  out <- matrix(NA, ncol = n_resp_patterns, nrow = n_resp_patterns)
  on_cov <- 1
  for (l in 1:length(resp_pattern_list)) {
    n_levels <- nrow(num_levels_p_cov[l])
    # TODO das erzeugt noch nicht das richtige muster hier, hier noch mal dran gehen und gucken
    # dass ich die matrix wie auf meinem zettel erzeuge
    if (l == 1) {
      out[1:n_levels, 1:n_levels] <- resp_pattern_list[[l]]
    } else {
      pos <- sum(num_levels_p_cov[1:l-1])+1:n_levels
      out[pos, pos] <- resp_pattern_list[[l]]
        
    }
    on_cov <- on_cov + 1
  }
}