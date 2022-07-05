# ridge_penalty ------------------------------------------------------------------------

ridge_penalty <- function(alphas, lambda) {
  lambda * sum(alphas^2)
}

# lasso_penalty ------------------------------------------------------------------------

lasso_penalty <- function(alphas, lambda) {
  lambda * sum(abs(alphas))
}

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

# parse_model ------------------------------------------------------------------------

parse_model <- function(model, data, data_long, person_id) {
  
  model_clean <- gsub("\n", "", model)
  model_parts <- unlist(strsplit(model_clean, "[;]"))
  # item data
  model_theta <- model_parts[grepl("theta =~", model_parts) | grepl("theta=~",model_parts)]
  item_names <- gsub("theta", "", model_theta)
  item_names <- gsub("=~", "", item_names)
  item_names <- trimws(unlist(strsplit(item_names, "[+]")))
  if (data_long) {
    item_names_long <- unique(unlist(strsplit(item_names, "[(]")))
    responses <- item_names_long[!grepl("::", item_names_long)]
    item_names_long <- item_names_long[grepl("::", item_names_long)]
    item_names_long <- unique(unlist(strsplit(item_names_long, "[::]")))
    item_names_long <- item_names_long[!("" == item_names_long)]
    item_id <- item_names_long[!grepl(")", item_names_long)]
    item_names <- gsub(")", "", item_names_long[grepl(")", item_names_long)])
    data_wide <- pivot_wider(
      data,
      id_cols = !!enquo(person_id),
      names_from = !!enquo(item_id),
      values_from = !!enquo(responses)
    )
    item_data <- data_wide[,item_names]
  } else {
    item_data <- data[,item_names, drop = FALSE]
  }
  if (length(model_parts) > 1) {
    # constraints
    # fixing parameters
    fixed_params <- any(grepl("=:", model_parts))
    if (fixed_params) {
      model_fixed <- model_parts[grepl("=:", model_parts)]
      # atm i can only fix log_nus and alphas
      if (any(grepl("alpha", model_parts))) {
        fixed_alphas <- model_fixed[grepl("alpha", model_fixed)]
        fixed_alphas <- gsub("alpha[0-9]+=:", "", fixed_alphas)
        fixed_alphas <- gsub("alpha[0-9]+ =:", "", fixed_alphas)
        fixed_alphas <- as.numeric(trimws(fixed_alphas))
      } else {
        fixed_alphas <- NULL
      }
      if (any(grepl("log_nu", model_parts))) {
        fixed_log_disps <- model_fixed[grepl("log_nu", model_fixed)]
        fixed_log_disps <- gsub("log_nu[0-9]+=:", "", fixed_log_disps)
        fixed_log_disps <- gsub("log_nu[0-9]+ =:", "", fixed_log_disps)
        fixed_log_disps <- as.numeric(trimws(fixed_log_disps))
      } else {
        fixed_log_disps <- NULL
      }
    } else {
      fixed_alphas <- NULL
      fixed_log_disps <- NULL
    }
    # equal parameters
    equal_params <- any((grepl("~ 1", model_parts) | grepl("~1", model_parts)) &
                          !grepl("[+]", model_parts))
    if (equal_params) {
      model_equal <- model_parts[(grepl("~ 1", model_parts) | 
                                    grepl("~1", model_parts)) &
                                   !grepl("[+]", model_parts)]
      # # atm i can only constrain log_nus and alphas to be equal across items
      if (any(grepl("alphas", model_equal))) {
        equal_alphas <- TRUE
      } else {
        equal_alphas <- FALSE
      }
      if (any(grepl("log_nus", model_equal))) {
        equal_log_disps <- TRUE
      } else {
        equal_log_disps <- FALSE
      }
    } else {
      equal_alphas <- FALSE
      equal_log_disps <- FALSE
    }
    # covariates
    # person covariates
    model_pcov <- model_parts[grepl("thetas ~", model_parts) | grepl("thetas~", model_parts)]
    if (length(model_pcov) > 0) {
      p_cov <- gsub("thetas", "", model_pcov)
      p_cov <- gsub("~", "", p_cov)
      p_cov <- trimws(unlist(strsplit(p_cov, "[+]")))
      p_cov <- p_cov[!p_cov == "1"]
      # TODO so far this will only work with data in wide format, so maybe add
      # possibility of long format here
      p_covariates <- data[, p_cov, drop = FALSE]
      if (all(sapply(p_covariates, function(x){is(x, "factor")}))) {
        # atm, we can only do p_cov_cat = TRUE when all p_covariates are factors
        p_cov_cat <- TRUE
        p_covariates <- lapply(
          p_covariates, 
          function(x){model.matrix(~x)[,-1,drop=FALSE]}
        )
        # if (length(p_cov) == 1) {
        #   p_covariates <- list(p_covariates)
        #   names(p_covariates) <- p_cov
        # }
        p_cov_names <- names(p_covariates)
        p_cov_levels <- unlist(lapply(p_covariates, ncol))+1
        p_covariates <- do.call(cbind, p_covariates)
        p_cov_names <- rep(p_cov_names, p_cov_levels-1)
        p_cov_levels_long <- unlist(sapply(p_cov_levels, function(x){2:x}))
        colnames(p_covariates) <- paste0(p_cov_names, ".", p_cov_levels_long)
      } else {
        p_cov_cat <- FALSE
        p_cov_levels <- NULL
      }
    } else {
      p_covariates <- NULL
      p_cov_cat <- FALSE
      p_cov_levels <- NULL
    }
    # item covariates
    model_icov <- model_parts[
      ((grepl("alphas ~", model_parts) | grepl("alphas~", model_parts)) |
        (grepl("deltas ~", model_parts) | grepl("deltas~", model_parts)) |
        grepl("log_nus ~", model_parts) | grepl("log_nus~", model_parts)) &
        !grepl("~ *1 *$", model_parts)
    ]
    if (length(model_icov) > 0) {
      i_cov_on <- trimws(gsub("~.*", "", model_icov))
      i_cov_on <- gsub("alphas", "alpha", i_cov_on)
      i_cov_on <- gsub("deltas", "delta", i_cov_on)
      i_cov_on <- gsub("log_nus", "log_disp", i_cov_on)
      # atm we can only have the same covariates on all item parameters
      i_cov <- model_icov[1]
      i_cov <- gsub("alphas", "", i_cov)
      i_cov <- gsub("deltas", "", i_cov)
      i_cov <- gsub("log_nus", "", i_cov)
      i_cov <- gsub("~", "", i_cov)
      i_cov <- trimws(unlist(strsplit(i_cov, "[+]")))
      i_cov <- i_cov[!i_cov == "1"]
      # for item covariates I am currently expecting long format
      # TODO implement a way for wide format here as well
      if (data_long) {
        # a column for each covariate with as many rows as we have items
        i_covariates <- unique(data[,c(i_cov, item_id), drop = FALSE])
        i_covariates <- i_covariates[,colnames(i_covariates) != item_id, drop = FALSE]
        rownames(i_covariates) <- 1:nrow(i_covariates)
        # TODO handling einfuegen fuer kategoriale kovaraiten und wie ich dafuer dann
        # hier dummy variablen erstelle
      } 
    } else {
      i_covariates <- NULL
      i_cov_on <- NULL
    }
  } else { # just a 2pcmpm with no constraints and no covariates
    fixed_alphas <- NULL
    fixed_log_disps <- NULL
    equal_alphas <- FALSE
    equal_log_disps <- FALSE
    p_covariates <- NULL
    i_covariates <- NULL
    i_cov_on <- NULL
    p_cov_cat <- FALSE
    p_cov_levels <- NULL
  }
  # otherwise we only want to fit a full 2pcmp with no covariates
  
  # prepare output
  out <- list(
    item_data = item_data,
    fixed_alphas = fixed_alphas,
    fixed_log_disps = fixed_log_disps,
    equal_alphas = equal_alphas,
    equal_log_disps = equal_log_disps,
    p_covariates = p_covariates,
    i_covariates = i_covariates,
    i_cov_on = i_cov_on,
    p_cov_cat = p_cov_cat,
    p_cov_levels = p_cov_levels
  )
  
  return(out)
}




