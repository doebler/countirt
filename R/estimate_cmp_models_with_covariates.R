
# estep_cmp_with_cov ---------------------------------------------------------------------
estep_cmp_with_cov <- function(data, item_params, 
                               p_covariates, i_covariates,
                               weights_and_nodes,
                               i_cov_on = c("alpha", "delta", "log_disp"),
                               p_cov_cat = TRUE,
                               resp_patterns_matrix = NULL) {
  # p_covariates is a matrix with the person covariates (MxP)
  # i_covariates is a matrix with the item covariates (M+I)
  
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # alphas is a scalar if we have covariates on alpha
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # in case of item covariates on delta, deltas will be ust a scalar
  log_disps <- item_params[grepl("log_disp", names(item_params)) &
                             !grepl("beta", names(item_params))]
  disps <- exp(log_disps)
  # in case of item covariates on nu, disps will be a scalar
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      PPs <- estep_cmp_with_pcov_cat_cpp(
        data = as.matrix(data),
        alphas = alphas,
        deltas = deltas,
        disps = disps,
        betas = betas_p,
        p_cov_data = as.matrix(p_covariates),
        resp_pattern = resp_patterns_matrix,
        nodes = weights_and_nodes$x,
        weights = weights_and_nodes$w,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_logZ_long = grid_logZ_long,
        grid_log_lambda_long = grid_log_lambda_long,
        max_mu = 200,
        min_mu = 0.001
      )
      # note: within the estep function, we work with response patterns here but
      # we then already select the correct response pattern for each person, so 
      # as always, we return a NxK matrix here
    } else {
      PPs <- estep_cmp_with_pcov_cpp(
        data = as.matrix(data),
        alphas = alphas,
        deltas = deltas,
        disps = disps,
        betas = betas_p,
        p_cov_data = as.matrix(p_covariates),
        nodes = weights_and_nodes$x,
        weights = weights_and_nodes$w,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_logZ_long = grid_logZ_long,
        grid_log_lambda_long = grid_log_lambda_long,
        max_mu = 200,
        min_mu = 0.001
      )
    }
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have the covaraites
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        PPs <- estep_cmp_with_icov_delta_cpp(
          data = as.matrix(data),
          alphas = alphas,
          delta = deltas,
          disps = disps,
          betas = betas_i,
          i_cov_data = as.matrix(i_covariates),
          nodes = weights_and_nodes$x,
          weights = weights_and_nodes$w,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_logZ_long = grid_logZ_long,
          grid_log_lambda_long = grid_log_lambda_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "alpha") {
        PPs <- estep_cmp_with_icov_alpha_cpp(
          data = as.matrix(data),
          alpha = alphas,
          deltas = deltas,
          disps = disps,
          betas = betas_i,
          i_cov_data = as.matrix(i_covariates),
          nodes = weights_and_nodes$x,
          weights = weights_and_nodes$w,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_logZ_long = grid_logZ_long,
          grid_log_lambda_long = grid_log_lambda_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "log_disp") {
        PPs <- estep_cmp_with_icov_nu_cpp(
          data = as.matrix(data),
          alphas = alphas,
          deltas = deltas,
          disp = disps,
          betas = betas_i,
          i_cov_data = as.matrix(i_covariates),
          nodes = weights_and_nodes$x,
          weights = weights_and_nodes$w,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_logZ_long = grid_logZ_long,
          grid_log_lambda_long = grid_log_lambda_long,
          max_mu = 200,
          min_mu = 0.001,
          max_nu = 50,
          min_nu = 0.001
        )
      }
    } else if (length(i_cov_on) == 3) {
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
      
      PPs <- estep_cmp_with_icov_all_cpp(
        data = as.matrix(data),
        alpha = alphas,
        delta = deltas,
        disp = disps,
        betas_alpha = betas_i_alpha,
        betas_delta = betas_i_delta,
        betas_logdisp = betas_i_logdisp,
        i_cov_data = as.matrix(i_covariates),
        nodes = weights_and_nodes$x,
        weights = weights_and_nodes$w,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_logZ_long = grid_logZ_long,
        grid_log_lambda_long = grid_log_lambda_long,
        max_mu = 200,
        min_mu = 0.001,
        max_nu = 50,
        min_nu = 0.001
      )
    } else if (length(i_cov_on) == 2) {
      if ("log_disp" %in% i_cov_on) {
        if ("alpha" %in% i_cov_on) {
          # we have covariates on alpha and log_disp
          betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
          betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
          
          PPs <- estep_cmp_with_icov_alpha_nu_cpp(
            data = as.matrix(data),
            alpha = alphas,
            deltas = deltas,
            disp = disps,
            betas_alpha = betas_i_alpha,
            betas_logdisp = betas_i_logdisp,
            i_cov_data = as.matrix(i_covariates),
            nodes = weights_and_nodes$x,
            weights = weights_and_nodes$w,
            grid_mus = grid_mus,
            grid_nus = grid_nus,
            grid_logZ_long = grid_logZ_long,
            grid_log_lambda_long = grid_log_lambda_long,
            max_mu = 200,
            min_mu = 0.001,
            max_nu = 50,
            min_nu = 0.001
          )
        } else {
          # we have covariates on delta and log_disp
          betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
          betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
          
          PPs <- estep_cmp_with_icov_delta_nu_cpp(
            data = as.matrix(data),
            alphas = alphas,
            delta = deltas,
            disp = disps,
            betas_delta = betas_i_delta,
            betas_logdisp = betas_i_logdisp,
            i_cov_data = as.matrix(i_covariates),
            nodes = weights_and_nodes$x,
            weights = weights_and_nodes$w,
            grid_mus = grid_mus,
            grid_nus = grid_nus,
            grid_logZ_long = grid_logZ_long,
            grid_log_lambda_long = grid_log_lambda_long,
            max_mu = 200,
            min_mu = 0.001,
            max_nu = 50,
            min_nu = 0.001
          )
        }
      } else {
        # we have covariates on alpha and delta
        betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        
        PPs <- estep_cmp_with_icov_alpha_delta_cpp(
          data = as.matrix(data),
          alpha = alphas,
          delta = deltas,
          disps = disps,
          betas_alpha = betas_i_alpha,
          betas_delta = betas_i_delta,
          i_cov_data = as.matrix(i_covariates),
          nodes = weights_and_nodes$x,
          weights = weights_and_nodes$w,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_logZ_long = grid_logZ_long,
          grid_log_lambda_long = grid_log_lambda_long,
          max_mu = 200,
          min_mu = 0.001
        )
      }
    }
  }
  
  return(PPs)
}

# grad_cmp_with_cov ----------------------------------------------------------------------
grad_cmp_with_cov <- function(item_params, PPs, weights_and_nodes, data, 
                              p_covariates, i_covariates,
                              i_cov_on = c("alpha", "delta", "log_disp"),
                              p_cov_cat = TRUE,
                              resp_patterns_matrix = NULL) {
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # alphas is a scalar if we have item covariates on alpha
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # note that for item covariates on delta, deltas is just a scalar
  log_disps <- item_params[grepl("log_disp", names(item_params)) &
                             !grepl("beta", names(item_params))]
  disps <- exp(log_disps)
  # disps is a scalar if we have item covariates on nu
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      grads <- grad_cmp_with_pcov_cat_cpp(
        alphas = alphas,
        deltas = deltas,
        disps = disps,
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        resp_pattern = resp_patterns_matrix,
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    } else {
      grads <- grad_cmp_with_pcov_cpp(
        alphas = alphas,
        deltas = deltas,
        disps = disps,
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    }
  } else if (!is.null(i_covariates)) { 
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_cpp(
          alphas = alphas,
          delta = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      } else if (i_cov_on == "alpha") {
        grads <- grad_cmp_with_icov_alpha_cpp(
          alpha = alphas,
          deltas = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      } else if (i_cov_on == "log_disp") {
        grads <- grad_cmp_with_icov_nu_cpp(
          alphas = alphas,
          deltas = deltas,
          disp = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001,
          max_nu = 50,
          min_nu = 0.001)
      }
    } else if (length(i_cov_on) == 3) {
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
      
      grads <- grad_cmp_with_icov_all_cpp(
        alpha = alphas,
        delta = deltas,
        disp = disps,
        betas_alpha = betas_i_alpha,
        betas_delta = betas_i_delta,
        betas_logdisp = betas_i_logdisp,
        data = as.matrix(data),
        i_cov_data = as.matrix(i_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001,
        max_nu = 50,
        min_nu = 0.001)
    } else if (length(i_cov_on) == 2) {
      if ("log_disp" %in% i_cov_on) {
        if ("alpha" %in% i_cov_on) {
          # we have covariates on log_disp and alpha together
          betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
          betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
          
          grads <- grad_cmp_with_icov_alpha_nu_cpp(
            alpha = alphas,
            deltas = deltas,
            disp = disps,
            betas_alpha = betas_i_alpha,
            betas_logdisp = betas_i_logdisp,
            data = as.matrix(data),
            i_cov_data = as.matrix(i_covariates),
            PPs = PPs,
            nodes = weights_and_nodes$x,
            grid_mus = grid_mus,
            grid_nus = grid_nus,
            grid_cmp_var_long = grid_cmp_var_long,
            grid_log_lambda_long = grid_log_lambda_long,
            grid_logZ_long = grid_logZ_long,
            max_mu = 200,
            min_mu = 0.001,
            max_nu = 50,
            min_nu = 0.001)
        } else {
          # we have covariates on log_disp and delta together
          betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
          betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
          
          grads <- grad_cmp_with_icov_delta_nu_cpp(
            alphas = alphas,
            delta = deltas,
            disp = disps,
            betas_delta = betas_i_delta,
            betas_logdisp = betas_i_logdisp,
            data = as.matrix(data),
            i_cov_data = as.matrix(i_covariates),
            PPs = PPs,
            nodes = weights_and_nodes$x,
            grid_mus = grid_mus,
            grid_nus = grid_nus,
            grid_cmp_var_long = grid_cmp_var_long,
            grid_log_lambda_long = grid_log_lambda_long,
            grid_logZ_long = grid_logZ_long,
            max_mu = 200,
            min_mu = 0.001,
            max_nu = 50,
            min_nu = 0.001)
        }
      } else {
        # we have covariates on alpha and delta together
        betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        
        grads <- grad_cmp_with_icov_alpha_delta_cpp(
          alpha = alphas,
          delta = deltas,
          disps = disps, 
          betas_alpha = betas_i_alpha,
          betas_delta = betas_i_delta,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      }
    }
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  }
  
  return(grads)
}

# grad_cmp_with_cov_fixdisps----------------------------------------------------------
grad_cmp_with_cov_fixdisps <- function(item_params, PPs, weights_and_nodes, 
                                       data, fix_disps,
                                       p_covariates, i_covariates,
                                       i_cov_on = c("alpha", "delta"),
                                       p_cov_cat = TRUE,
                                       resp_patterns_matrix = NULL) {

  alphas <- item_params[grepl("alpha", names(item_params))]
  # note that alphas is a scalar if we have item covariates on alpha
  deltas <- item_params[grepl("delta", names(item_params))]
  # note that deltas is a scalar if we have item covariates on delta
  disps <- fix_disps
  n_items <- ncol(data)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      grads <- grad_cmp_with_pcov_cat_fixdisps_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        resp_pattern = resp_patterns_matrix,
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    } else {
      grads <- grad_cmp_with_pcov_fixdisps_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    }
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameters we have the covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_fixdisps_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "alpha") {
        grads <- grad_cmp_with_icov_alpha_fixdisps_cpp(
          alpha = alphas, 
          deltas = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } 
      # note: we can't include item covariates on log nu if we fix dispersions
      # to a specific value
    } else if (length(i_cov_on) == 2) {
      # if we have the constraint of fixed disps, we can only have item covariates on 
      # alpha and delta
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      
      grads <- grad_cmp_with_icov_alpha_delta_fixdisps_cpp(
        alpha = alphas,
        delta = deltas,
        disps = disps,
        betas_alpha = betas_i_alpha,
        betas_delta = betas_i_delta,
        data = as.matrix(data),
        i_cov_data = as.matrix(i_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    }
    # note that we can't have item covariates on all item parameters if we have the
    # constraint of fixed disps as predicting all parameters, i.e., including log_disp,
    # through item covariates implies different disps for items with different
    # values on the item covariates
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  return(grads)
}

# grad_cmp_with_cov_fixalphas-------------------------------------------------------
grad_cmp_with_cov_fixalphas <- function(item_params, PPs, weights_and_nodes, 
                                        data, fix_alphas,
                                        p_covariates, i_covariates,
                                        i_cov_on = c("delta", "log_disp"),
                                        p_cov_cat = TRUE, 
                                        resp_patterns_matrix = NULL) {
  
  alphas <- fix_alphas
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # note that deltas is a scalar if we have item covariates on delta
  log_disps <- item_params[grepl("log_disp", names(item_params)) &
                             !grepl("beta", names(item_params))]
  disps <- exp(log_disps)
  # note that disps is a scalar if we have item covariates on nu
  n_items <- ncol(data)
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      grads <- grad_cmp_with_pcov_cat_fixalphas_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        resp_pattern = resp_patterns_matrix,
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    } else {
      grads <- grad_cmp_with_pcov_fixalphas_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    }
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have the covariates
    if (length(i_cov_on) == 1) {
      # note: if we have the constraint that alphas are fixed to specfic values
      # then we can't predict alphas through covariates
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_fixalphas_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "log_disp") {
        grads <- grad_cmp_with_icov_nu_fixalphas_cpp(
          alphas = alphas, 
          deltas = deltas, 
          disp = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001,
          max_nu = 50,
          min_nu = 0.001
        )
      }
    } else if (length(i_cov_on) == 2) {
      # with the constraint same alphas, we can only have covariates on
      # delta and log disp at the same time
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
      
      grads <- grad_cmp_with_icov_delta_nu_fixalphas_cpp(
        alphas = alphas,
        delta = deltas,
        disp = disps,
        betas_delta = betas_i_delta,
        betas_logdisp = betas_i_logdisp,
        data = as.matrix(data),
        i_cov_data = as.matrix(i_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001,
        max_nu = 50,
        min_nu = 0.001)
    }
    # note that we can't have item covariates on all item parameters if we have the
    # constraint of same alphas as predicting all parameters, i.e., including alpha,
    # through item covariates implies different alphas for items with different
    # values on the item covariates
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  return(grads)
}

# grad_cmp_with_cov_samedisps ---------------------------------------------------------
grad_cmp_with_cov_samedisps <- function(item_params, PPs, 
                                        weights_and_nodes, data,
                                        p_covariates, i_covariates,
                                        i_cov_on = c("alpha", "delta"),
                                        p_cov_cat = TRUE,
                                        resp_patterns_matrix = NULL) {
  
  alphas <- item_params[grepl("alpha", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # note that alphas is a scalar for item covariates on alpha
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # note that deltas is a scalar for item covariates on delta
  n_items <- ncol(data)
  log_disp <- item_params[grepl("log_disp", names(item_params))]
  disps <- exp(rep(log_disp, n_items))
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      grads <- grad_cmp_with_pcov_cat_samedisps_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        resp_pattern = resp_patterns_matrix,
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    } else {
      grads <- grad_cmp_with_pcov_samedisps_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    }
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_samedisps_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(c_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "alpha") {
        grads <- grad_cmp_with_icov_alpha_samedisps_cpp(
          alpha = alphas, 
          deltas = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(c_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      }
      # note: we can't have covariates on log_nu if we have the constraint of
      # same disps as covaraites would have different values for the different
      # items and would then imply different disps
    } else if (length(i_cov_on) == 2) {
      # if we have the constraint of same disps, we can only have item covariates on 
      # alpha and delta
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      
      grads <- grad_cmp_with_icov_alpha_delta_samedisps_cpp(
        alpha = alphas,
        delta = deltas,
        disps = disps,
        betas_alpha = betas_i_alpha,
        betas_delta = betas_i_delta,
        data = as.matrix(data),
        i_cov_data = as.matrix(i_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    }
    # note that we can't have item covariates on all item parameters if we have the
    # constraint of same disps as predicting all parameters, i.e., including log_disp,
    # through item covariates implies different disps for items with different
    # values on the item covariates
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 

  
  return(grads)
}

# grad_cmp_with_cov_samealphas -----------------------------------------------------
grad_cmp_with_cov_samealphas <- function(item_params, PPs, 
                                         weights_and_nodes, data,
                                         p_covariates, i_covariates,
                                         i_cov_on = c("delta", "log_disp"),
                                         p_cov_cat = TRUE, 
                                         resp_patterns_matrix = NULL) {
  
  deltas <- item_params[grepl("delta", names(item_params)) & 
                          !grepl("beta", names(item_params))]
  # note that for item covariates on delta, deltas is a scalar
  n_items <- ncol(data)
  alpha <- item_params[grepl("alpha", names(item_params))]
  alphas <- rep(alpha, n_items)
  # note: we can't have covariates on alpha if we have the constraint same_alpha
  # as we would have different values for the items on the different covariates,
  # implying different alphas
  log_disps <- item_params[grepl("log_disp", names(item_params)) &
                             !grepl("beta", names(item_params))]
  disps <- exp(log_disps)
  # note that for item covariates on log disp, disps is a scalar
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      grads <- grad_cmp_with_pcov_cat_samealphas_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        resp_pattern = resp_patterns_matrix,
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    } else {
      grads <- grad_cmp_with_pcov_samealphas_cpp(
        alphas = alphas, 
        deltas = deltas, 
        disps = disps, 
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001
      )
    }
  } else if (!is.null(i_covariates)) {
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on) == 1) {
      # note: we can't have covariates on alpha if we have the constraint same_alpha
      # as we would have different values for the items on the different covariates,
      # implying different alphas
      if (i_cov_on == "delta") {
        grads <- grad_cmp_with_icov_delta_samealphas_cpp(
          alphas = alphas, 
          delta = deltas, 
          disps = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001
        )
      } else if (i_cov_on == "log_disp") {
        grads <- grad_cmp_with_icov_nu_samealphas_cpp(
          alphas = alphas, 
          deltas = deltas, 
          disp = disps, 
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001,
          max_nu = 50,
          min_nu = 0.001
        )
      }
    } else if (length(i_cov_on) == 2) {
      # with the constraint same alphas, we can only have covariates on
      # delta and log disp at the same time
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
      
      grads <- grad_cmp_with_icov_delta_nu_samealphas_cpp(
        alphas = alphas,
        delta = deltas,
        disp = disps,
        betas_delta = betas_i_delta,
        betas_logdisp = betas_i_logdisp,
        data = as.matrix(data),
        i_cov_data = as.matrix(i_covariates),
        PPs = PPs,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001,
        max_nu = 50,
        min_nu = 0.001)
    }
    # note that we can't have item covariates on all item parameters if we have the
    # constraint of same alphas as predicting all parameters, i.e., including alpha,
    # through item covariates implies different alphas for items with different
    # values on the item covariates
  }
  
  if (any(is.na(grads))) {
    stop("Gradient contained NA", paste0(grads, collapse = ","),
         paste0(item_params, collapse = ","))
  } 
  
  return(grads)
}

# ell_cmp_with_cov -------------------------------------------------------------------
ell_cmp_with_cov <- function(item_params, PPs, weights_and_nodes, 
                             data, p_covariates, i_covariates,
                             i_cov_on = c("alpha", "delta", "log_disp"),
                             p_cov_cat = TRUE,
                             resp_patterns_matrix = NULL) {
  # ell without any restrictions
  
  # prep item parameters
  alphas <- item_params[grepl("alpha", names(item_params)) & 
                          !grepl("beta", names(item_params))]
  # note that akphas is a scalar if we have covariates on alpha
  deltas <- item_params[grepl("delta", names(item_params)) &
                          !grepl("beta", names(item_params))]
  # note that deltas is a scalar if we have item covariates on delta
  log_disps <- item_params[grepl("log_disp", names(item_params)) 
                           & !grepl("beta", names(item_params))]
  disps <- exp(log_disps)
  # disps will be a scalar if we have covariates on log disps
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      ell <- ell_cmp_with_pcov_cat_cpp(
        alphas = alphas,
        deltas = deltas,
        disps = disps,
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        resp_pattern = resp_patterns_matrix,
        PPs = PPs,
        weights = weights_and_nodes$w,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    } else {
      ell <- ell_cmp_with_pcov_cpp(
        alphas = alphas,
        deltas = deltas,
        disps = disps,
        betas = betas_p,
        data = as.matrix(data),
        p_cov_data = as.matrix(p_covariates),
        PPs = PPs,
        weights = weights_and_nodes$w,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001)
    }
  } else if (!is.null(i_covariates)) { 
    # distinguish between on which item parameter we have covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        ell <- ell_cmp_with_icov_delta_cpp(
          alphas = alphas,
          delta = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          weights = weights_and_nodes$w,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      } else if (i_cov_on == "alpha") {
        ell <- ell_cmp_with_icov_alpha_cpp(
          alpha = alphas,
          deltas = deltas,
          disps = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          weights = weights_and_nodes$w,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      } else if (i_cov_on == "log_disp") {
        ell <- ell_cmp_with_icov_nu_cpp(
          alphas = alphas,
          deltas = deltas,
          disp = disps,
          betas = betas_i,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          weights = weights_and_nodes$w,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001,
          max_nu = 50,
          min_nu = 0.001)
      }
    } else if (length(i_cov_on) == 3) {
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
      
      ell <- ell_cmp_with_icov_all_cpp(
        alpha = alphas,
        delta = deltas,
        disp = disps,
        betas_alpha = betas_i_alpha,
        betas_delta = betas_i_delta,
        betas_logdisp = betas_i_logdisp,
        data = as.matrix(data),
        i_cov_data = as.matrix(i_covariates),
        PPs = PPs,
        weights = weights_and_nodes$w,
        nodes = weights_and_nodes$x,
        grid_mus = grid_mus,
        grid_nus = grid_nus,
        grid_cmp_var_long = grid_cmp_var_long,
        grid_log_lambda_long = grid_log_lambda_long,
        grid_logZ_long = grid_logZ_long,
        max_mu = 200,
        min_mu = 0.001,
        max_nu = 50,
        min_nu = 0.001)
    } else if (length(i_cov_on) == 2) {
      if ("log_disp" %in% i_cov_on) {
        if ("alpha" %in% i_cov_on) {
          # covariates on alpha and log disp
          betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
          betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
          
          ell <- ell_cmp_with_icov_alpha_nu_cpp(
            alpha = alphas,
            deltas = deltas,
            disp = disps,
            betas_alpha = betas_i_alpha,
            betas_logdisp = betas_i_logdisp,
            data = as.matrix(data),
            i_cov_data = as.matrix(i_covariates),
            PPs = PPs,
            weights = weights_and_nodes$w,
            nodes = weights_and_nodes$x,
            grid_mus = grid_mus,
            grid_nus = grid_nus,
            grid_cmp_var_long = grid_cmp_var_long,
            grid_log_lambda_long = grid_log_lambda_long,
            grid_logZ_long = grid_logZ_long,
            max_mu = 200,
            min_mu = 0.001,
            max_nu = 50,
            min_nu = 0.001)
        } else {
          # covariates on delta and log disp
          betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
          betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
          
          ell <- ell_cmp_with_icov_delta_nu_cpp(
            alphas = alphas,
            delta = deltas,
            disp = disps,
            betas_delta = betas_i_delta,
            betas_logdisp = betas_i_logdisp,
            data = as.matrix(data),
            i_cov_data = as.matrix(i_covariates),
            PPs = PPs,
            weights = weights_and_nodes$w,
            nodes = weights_and_nodes$x,
            grid_mus = grid_mus,
            grid_nus = grid_nus,
            grid_cmp_var_long = grid_cmp_var_long,
            grid_log_lambda_long = grid_log_lambda_long,
            grid_logZ_long = grid_logZ_long,
            max_mu = 200,
            min_mu = 0.001,
            max_nu = 50,
            min_nu = 0.001)
        }
      } else {
        # covariates on alpha and delta together
        betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        
        ell <- ell_cmp_with_icov_all_cpp(
          alpha = alphas,
          delta = deltas,
          disps = disps,
          betas_alpha = betas_i_alpha,
          betas_delta = betas_i_delta,
          data = as.matrix(data),
          i_cov_data = as.matrix(i_covariates),
          PPs = PPs,
          weights = weights_and_nodes$w,
          nodes = weights_and_nodes$x,
          grid_mus = grid_mus,
          grid_nus = grid_nus,
          grid_cmp_var_long = grid_cmp_var_long,
          grid_log_lambda_long = grid_log_lambda_long,
          grid_logZ_long = grid_logZ_long,
          max_mu = 200,
          min_mu = 0.001)
      }
    }
  }
  
  return(ell)
}

# em_cycle_cmp_with_cov ---------------------------------------------------------------------
em_cycle_cmp_with_cov <- function(data, item_params, weights_and_nodes,
                                  p_covariates, i_covariates,
                                  i_cov_on = c("alpha", "delta", "log_disp"),
                                  p_cov_cat = TRUE,
                                  resp_patterns_matrix = NULL,
                                  ctol_maxstep = 1e-8, m_method = "nleqslv",
                                  fix_disps = NULL, fix_alphas = NULL,
                                  same_disps = FALSE, same_alphas = FALSE) {
  
  if (is.null(fix_disps) & is.null(fix_alphas)) {
    if (!same_disps & !same_alphas) {
      # e step
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!same_disps & same_alphas) { 
      
      # prep for e step
      alpha <- item_params[grepl("alpha", names(item_params)) &
                             !grepl("beta", names(item_params))]
      deltas <- item_params[grepl("delta", names(item_params)) &
                              !grepl("beta", names(item_params))]
      n_items <- ncol(data)
      log_disps <- item_params[grepl("log_disp", names(item_params)) &
                                 !grepl("beta", names(item_params))]
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      
      # if we have the constraint of same alphas, we can either have covariates
      # on just delta or just log_disp, in which case we just have betas_i,
      # or we could have covariates on delta and log_disp together (or person covariates)
      # with the constraint, we can't have item covariates on all three parameters
      if (length(i_cov_on) == 2) {
        # must be covariates on delta and log_disp together
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
        item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps,  
                                  betas_i_delta, betas_i_logdisp)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params)) & 
                              !grepl("beta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params)) & 
                              !grepl("beta", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))]),
          names(item_params[grepl("beta_i_log_disp", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samealph <- c(rep(alpha, n_items), deltas, log_disps, betas_p, betas_i)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params))]),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }

      # e step
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_samealph,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_samealphas,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (same_disps & !same_alphas) {

      # prep the parameters for the e-step
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
        item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items), 
                                  betas_i_alpha, betas_i_delta)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_i_alpha", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samedisp <- c(alphas, deltas, rep(log_disp, n_items), betas_p, betas_i)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }
      
      # e step
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_samedisp,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_samedisps,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  } else {
    if (!is.null(fix_disps)) {
      
      # prep for e step
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
        item_params_samedisp <- c(alphas, deltas, log(fix_disps), 
                                  betas_i_alpha, betas_i_delta)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_i_alpha", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samedisp <- c(alphas, deltas, log(fix_disps), betas_p, betas_i)
        names(item_params_samedisp) <- c(
          names(item_params[grepl("alpha", names(item_params))]),
          names(item_params[grepl("delta", names(item_params))]),
          paste0("log_disp", 1:n_items),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }
      
      # e step
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_fixdisps,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_fixdisps,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_disps = fix_disps,
        control = list(xtol = ctol_maxstep)
      )$x
    } else if (!is.null(fix_alphas)) {
      
      # prep for e step
      deltas <- item_params[grepl("delta", names(item_params)) &
                              !grepl("beta", names(item_params))]
      n_items <- ncol(data)
      log_disps <- item_params[grepl("log_disp", names(item_params)) &
                                 !grepl("beta", names(item_params))]
      betas_p <- item_params[grepl("beta_p", names(item_params))]
      betas_i <- item_params[grepl("beta_i", names(item_params))]
      
      # if we have the constraint of fixed alphas, we can either have covariates
      # on just delta or just log_disp, in which case we just have betas_i,
      # or we could have covariates on delta and log_disp together (or person covariates)
      # but with the constraint, we can't have item covaraites on all three parameters
      if (length(i_cov_on) == 2) {
        # must be covariates on delta and log_disp together
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        betas_i_logdisp <- betas_i[grepl("log_disp", names(betas_i))]
        item_params_samealph <- c(fix_alphas, deltas, log_disps,  
                                  betas_i_delta, betas_i_logdisp)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params)) & 
                              !grepl("beta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params))]),
          names(item_params[grepl("beta_i_delta", names(item_params))]),
          names(item_params[grepl("beta_i_log_disp", names(item_params))])
        )
      } else {
        # either we have person covariates or just item covariates on one parameter, so just beta_i
        item_params_samealph <- c(fix_alphas, deltas, log_disps, betas_p, betas_i)
        names(item_params_samealph) <- c(
          paste0("alpha", 1:n_items),
          names(item_params[grepl("delta", names(item_params))]),
          names(item_params[grepl("log_disp", names(item_params))]),
          names(item_params[grepl("beta_p", names(item_params))]),
          names(item_params[grepl("beta_i", names(item_params))])
        )
      }

      # e step 
      PPs <- estep_cmp_with_cov(
        data = data, 
        item_params = item_params_fixalphas,
        weights_and_nodes = weights_and_nodes,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix
      )
      # m step
      new_item_params <- nleqslv(
        x = item_params,
        fn = grad_cmp_with_cov_fixalphas,
        PPs = PPs,
        weights_and_nodes = weights_and_nodes,
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_alphas = fix_alphas,
        control = list(xtol = ctol_maxstep)
      )$x
    }
  }
  
  return(new_item_params)
}


# marg_ll_cmp_with_cov --------------------------------------------------------------------------

marg_ll_cmp_with_cov <- function(data, item_params, weights_and_nodes, 
                                 p_covariates, i_covariates, 
                                 i_cov_on = c("alpha", "delta", "log_disp"),
                                 p_cov_cat = TRUE,
                                 resp_patterns_matrix = NULL,
                                 fix_disps = NULL, fix_alphas = NULL, 
                                 same_disps = FALSE, same_alphas = FALSE) {
  n_items <- ncol(data)
  n_persons <- nrow(data)
  deltas <- item_params[grepl("delta", names(item_params)) & 
                          !grepl("beta", names(item_params))]
  # note that deltas will be a scalar if we have item covariates on delta
  if (is.null(fix_alphas)) {
    if (same_alphas) {
      alpha <- item_params[grepl("alpha", names(item_params))]
      alphas <- rep(alpha, n_items)
    } else {
      alphas <- item_params[grepl("alpha", names(item_params)) &
                              !grepl("beta", names(item_params))]
      # alphas will be a scalar if we have item covariates on alpha
    }
  } else {
     alphas <- fix_alphas
  }
  if (is.null(fix_disps)) {
    if (same_disps) {
      log_disp <- item_params[grepl("log_disp", names(item_params))]
      disps <- rep(exp(log_disp), n_items)
    } else {
      log_disps <- item_params[grepl("log_disp", names(item_params)) &
                                 !grepl("beta", names(item_params))]
      disps <- exp(log_disps)
      # disps will be a scalar if we have item covariates on log_disp
    }
  } else {
    disps <- fix_disps
  }
  betas_p <- item_params[grepl("beta_p", names(item_params))]
  betas_i <- item_params[grepl("beta_i", names(item_params))]
  
  if (!is.null(p_covariates)) {
    if (p_cov_cat) {
      # assume that we have previously dummy coded our covariates, so we only have 0
      # and 1's in our p_cov matrix
      ll <- marg_ll_cmp_with_pcov_cat_cpp(data = as.matrix(data),
                                      alphas = alphas,
                                      deltas = deltas, 
                                      disps = disps, 
                                      betas = betas_p,
                                      p_cov_data = as.matrix(p_covariates),
                                      resp_pattern = resp_patterns_matrix,
                                      nodes = weights_and_nodes$x,
                                      weights = weights_and_nodes$w,
                                      grid_mus = grid_mus,  
                                      grid_nus = grid_nus, 
                                      grid_logZ_long = grid_logZ_long,
                                      grid_log_lambda_long = grid_log_lambda_long,
                                      max_mu = 150,
                                      min_mu = 0.001)
    } else {
      ll <- marg_ll_cmp_with_pcov_cpp(data = as.matrix(data),
                                      alphas = alphas,
                                      deltas = deltas, 
                                      disps = disps, 
                                      betas = betas_p,
                                      p_cov_data = as.matrix(p_covariates),
                                      nodes = weights_and_nodes$x,
                                      weights = weights_and_nodes$w,
                                      grid_mus = grid_mus,  
                                      grid_nus = grid_nus, 
                                      grid_logZ_long = grid_logZ_long,
                                      grid_log_lambda_long = grid_log_lambda_long,
                                      max_mu = 150,
                                      min_mu = 0.001)
    }
  } else if (!is.null(i_covariates)) {
    # distinguish between on which parameters we have item covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        ll <- marg_ll_cmp_with_icov_delta_cpp(data = as.matrix(data),
                                        alphas = alphas,
                                        delta = deltas, 
                                        disps = disps, 
                                        betas = betas_i,
                                        i_cov_data = as.matrix(i_covariates),
                                        nodes = weights_and_nodes$x,
                                        weights = weights_and_nodes$w,
                                        grid_mus = grid_mus,  
                                        grid_nus = grid_nus, 
                                        grid_logZ_long = grid_logZ_long,
                                        grid_log_lambda_long = grid_log_lambda_long,
                                        max_mu = 150,
                                        min_mu = 0.001)
      } else if (i_cov_on == "alpha") {
        ll <- marg_ll_cmp_with_icov_alpha_cpp(data = as.matrix(data),
                                              alpha = alphas,
                                              deltas = deltas, 
                                              disps = disps, 
                                              betas = betas_i,
                                              i_cov_data = as.matrix(i_covariates),
                                              nodes = weights_and_nodes$x,
                                              weights = weights_and_nodes$w,
                                              grid_mus = grid_mus,  
                                              grid_nus = grid_nus, 
                                              grid_logZ_long = grid_logZ_long,
                                              grid_log_lambda_long = grid_log_lambda_long,
                                              max_mu = 150,
                                              min_mu = 0.001)
      } else if (i_cov_on == "log_disp") {
        ll <- marg_ll_cmp_with_icov_nu_cpp(data = as.matrix(data),
                                              alphas = alphas,
                                              deltas = deltas, 
                                              disp = disps, 
                                              betas = betas_i,
                                              i_cov_data = as.matrix(i_covariates),
                                              nodes = weights_and_nodes$x,
                                              weights = weights_and_nodes$w,
                                              grid_mus = grid_mus,  
                                              grid_nus = grid_nus, 
                                              grid_logZ_long = grid_logZ_long,
                                              grid_log_lambda_long = grid_log_lambda_long,
                                              max_mu = 150,
                                              min_mu = 0.001,
                                           max_nu = 50,
                                           min_nu = 0.001)
      }
    } else if (length(i_cov_on) == 3) {
      betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
      betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
      betas_i_log_disp <- betas_i[grepl("log_disp", names(betas_i))]
      
      ll <- marg_ll_cmp_with_icov_all_cpp(data = as.matrix(data),
                                         alpha = alphas,
                                         delta = deltas, 
                                         disp = disps, 
                                         betas_alpha = betas_i_alpha,
                                         betas_delta = betas_i_delta,
                                         betas_logdisp = betas_i_log_disp,
                                         i_cov_data = as.matrix(i_covariates),
                                         nodes = weights_and_nodes$x,
                                         weights = weights_and_nodes$w,
                                         grid_mus = grid_mus,  
                                         grid_nus = grid_nus, 
                                         grid_logZ_long = grid_logZ_long,
                                         grid_log_lambda_long = grid_log_lambda_long,
                                         max_mu = 150,
                                         min_mu = 0.001,
                                         max_nu = 50,
                                         min_nu = 0.001)
    } else if (length(i_cov_on) == 2) {
      if ("log_disp" %in% i_cov_on) {
        betas_i_log_disp <- betas_i[grepl("log_disp", names(betas_i))]
        if ("alpha" %in% i_cov_on) {
          betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
          ll <- marg_ll_cmp_with_icov_alpha_nu_cpp(data = as.matrix(data),
                                                   alpha = alphas,
                                                   deltas = deltas, 
                                                   disp = disps, 
                                                   betas_alpha = betas_i_alpha,
                                                   betas_logdisp = betas_i_log_disp,
                                                   i_cov_data = as.matrix(i_covariates),
                                                   nodes = weights_and_nodes$x,
                                                   weights = weights_and_nodes$w,
                                                   grid_mus = grid_mus,  
                                                   grid_nus = grid_nus, 
                                                   grid_logZ_long = grid_logZ_long,
                                                   grid_log_lambda_long = grid_log_lambda_long,
                                                   max_mu = 150,
                                                   min_mu = 0.001,
                                                   max_nu = 50,
                                                   min_nu = 0.001)
        } else {
          # then we have covraitaes on delta and log_disp
          betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
          ll <- marg_ll_cmp_with_icov_delta_nu_cpp(data = as.matrix(data),
                                                   alphas = alphas,
                                                   delta = deltas, 
                                                   disp = disps, 
                                                   betas_delta = betas_i_delta,
                                                   betas_logdisp = betas_i_log_disp,
                                                   i_cov_data = as.matrix(i_covariates),
                                                   nodes = weights_and_nodes$x,
                                                   weights = weights_and_nodes$w,
                                                   grid_mus = grid_mus,  
                                                   grid_nus = grid_nus, 
                                                   grid_logZ_long = grid_logZ_long,
                                                   grid_log_lambda_long = grid_log_lambda_long,
                                                   max_mu = 150,
                                                   min_mu = 0.001,
                                                   max_nu = 50,
                                                   min_nu = 0.001)
        }
      } else {
        # then we have covariates on alpha and delta
        betas_i_alpha <- betas_i[grepl("alpha", names(betas_i))]
        betas_i_delta <- betas_i[grepl("delta", names(betas_i))]
        ll <- marg_ll_cmp_with_icov_alpha_delta_cpp(data = as.matrix(data),
                                                    alpha = alphas,
                                                    delta = deltas, 
                                                    disps = disps, 
                                                    betas_alpha = betas_i_alpha,
                                                    betas_delta = betas_i_delta,
                                                    i_cov_data = as.matrix(i_covariates),
                                                    nodes = weights_and_nodes$x,
                                                    weights = weights_and_nodes$w,
                                                    grid_mus = grid_mus, 
                                                    grid_nus = grid_nus, 
                                                    grid_logZ_long = grid_logZ_long,
                                                    grid_log_lambda_long = grid_log_lambda_long,
                                                    max_mu = 150,
                                                    min_mu = 0.001)
      }
    } 
  }
    
  return(ll)
}

# run_em_cmp_with_cov ----------------------------------------------------------------------
run_em_cmp_with_cov <- function(data, init_params, n_nodes, 
                                p_covariates, i_covariates, 
                                i_cov_on = c("alpha", "delta", "log_disp"),
                                p_cov_cat = TRUE,
                                num_levels_p_cov = NULL,
                                thres = Inf, prob = 0,
                                maxiter = 1000, convtol = 1e-5, ctol_maxstep = 1e-8,
                                m_method = "nleqslv", convcrit = "marglik",
                                fix_disps = NULL, fix_alphas = NULL,
                                same_disps = FALSE, same_alphas = FALSE) {

  # get nodes and weights for GH quadrature
  # weights_and_nodes <- gaussHermiteData(n_nodes)
  # weights_and_nodes$x <- weights_and_nodes$x * sqrt(2)
  # weights_and_nodes$w <- weights_and_nodes$w / sqrt(pi)
  weights_and_nodes<- quad_rule(n_nodes, thres = thres,prob = prob)

  new_params <- init_params
  conv <- FALSE
  iter <- 1

  new_ll <- 0
  marg_lls <- c()
  
  # prepare response patterns for categorical covariates
  if (!is.null(p_covariates) & p_cov_cat) {
    # create a possible response patterns matrix for the dummy coded covariates
    n_resp_patterns <- prod(num_levels_p_cov)
    # for each covariate, I first create a matrix with their possible response pattens
    # the first is always 0 erevrywhere and then from level2 to highest level resp.
    # one 1 and otherwise 0
    cov_list_resp_patterns <- lapply(num_levels_p_cov, get_resp_patterns_pcov_cat)
    resp_patterns_matrix <- make_resp_patterns_mat(
      cov_list_resp_patterns, n_resp_patterns, num_levels_p_cov
    )
  }

  print("Start estimation...")

  while (!isTRUE(conv) && (iter <= maxiter)) {
    print(paste0("Iteration: ", iter))
    old_params <- new_params
    new_params <- em_cycle_cmp_with_cov(
      data = data, 
      item_params = old_params, 
      weights_and_nodes = weights_and_nodes,
      p_covariates = p_covariates, 
      i_covariates = i_covariates,
      i_cov_on = i_cov_on,
      p_cov_cat = p_cov_cat,
      resp_patterns_matrix = resp_patterns_matrix,
      ctol_maxstep = ctol_maxstep,
      m_method = m_method,
      fix_disps = fix_disps, 
      fix_alphas = fix_alphas,
      same_disps = same_disps, 
      same_alphas = same_alphas
    )
    #print(new_params)

    # check for convergence
    if (convcrit == "marglik") {
      old_ll <- new_ll
      new_ll <-  marg_ll_cmp_with_cov(
        data = data,
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        p_covariates = p_covariates, 
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_disps = fix_disps, 
        fix_alphas = fix_alphas,
        same_disps = same_disps, 
        same_alphas = same_alphas)
      marg_lls[iter] <- new_ll
      #plot(marg_lls)
      #print(marg_lls)
      conv <- (abs(old_ll - new_ll) < convtol)
    } else {
      # convergence is to be assessed on parameter values, argument convcrit = "params"
      conv <- !any(abs(old_params - new_params) > convtol)
      marg_ll <- marg_ll_cmp_with_cov(
        data = data,
        item_params = new_params,
        weights_and_nodes = weights_and_nodes, 
        p_covariates = p_covariates, 
        i_covariates = i_covariates,
        i_cov_on = i_cov_on,
        p_cov_cat = p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        fix_disps = fix_disps, 
        fix_alphas = fix_alphas,
        same_disps = same_disps, 
        same_alphas = same_alphas)
      marg_lls[iter] <- marg_ll
      plot(marg_lls)
      print(marg_lls)
    }

    iter <- iter + 1
  }

  print("Done!")
  
  # model_vcov <- compute_vcov(
  #   item_params = new_params,
  #   weights_and_nodes = weights_and_nodes, 
  #   data = data
  # )
  # 
  # se_params <- se_from_vcov(model_vcov)

  out <- list(
    params = new_params,
#    se_params = se_params,
    iter = iter,
    conv = conv,
#    vcov = model_vcov,
    marg_ll = marg_lls
  )
  return(out)

}

# get_start_values_cmp_with_cov ---------------------------------------------------------------------

get_start_values_cmp_with_cov <- function(data, 
                                          p_covariates,
                                          i_covariates,
                                          nodes = 121, nsim = 1000,
                                          same_alpha = FALSE,
                                          i_cov_on = c("alpha", "delta", "log_disp")) {
  # p_covariates will work just like this, because poisson takes continuous or
  # categorical with dummy coding as it is and over than that we're not using it
  # here 
  
  # for CMP start values, we fit a Poisson model and get deltas and alphas 
  # and betas from there
  if (same_alpha) {
    if (!is.null(i_covariates)) {
      # just one alpha for all items
      # note that we can't have same alpha together with item covariates on alpha
      if (length(i_cov_on) == 1) {
        if (i_cov_on == "log_disp") {
          init_values_pois <- get_start_values_pois(
            data = data,
            same_alpha = TRUE
          )
          fit_pois <- run_em_poisson(
            data = data,
            init_params = init_values_pois,
            n_nodes = nodes,
            same_alpha = TRUE
          )
        } else {
          init_values_pois <- get_start_values_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            same_alpha = TRUE,
            i_cov_on = "delta"
          )
          fit_pois <- run_em_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            init_params = init_values_pois, 
            n_nodes = nodes,
            same_alpha = TRUE,
            i_cov_on = "delta"
          )
        }
      } else {
        # if I have the constrained of same alpha, i can't have covariates on all three
        # item parameters but only on two: log_disp and delta (also not on alpha and delta
        # or alpha and log_disp pairing)
        init_values_pois <- get_start_values_poisson_with_cov(
          data = data,
          p_covariates = p_covariates,
          i_covariates = i_covariates,
          same_alpha = TRUE,
          i_cov_on = "delta"
        )
        fit_pois <- run_em_poisson_with_cov(
          data = data,
          p_covariates = p_covariates,
          i_covariates = i_covariates,
          init_params = init_values_pois, 
          n_nodes = nodes,
          same_alpha = TRUE,
          i_cov_on = "delta"
        )
      } 
    } else if (!is.null(p_covariates)) {
      init_values_pois <- get_start_values_pois(
        data = data,
        same_alpha = TRUE
      )
      fit_pois <- run_em_poisson(
        data = data,
        init_params = init_values_pois,
        n_nodes = nodes,
        same_alpha = TRUE
      )
    }
  } else {
    if (!is.null(i_covariates)) {
      # no constraint on alpha
      if (length(i_cov_on) == 1) {
        if (i_cov_on == "log_disp") {
          init_values_pois <- get_start_values_pois(
            data = data
          )
          fit_pois <- run_em_poisson(
            data = data,
            init_params = init_values_pois,
            n_nodes = nodes
          )
        } else {
          init_values_pois <- get_start_values_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            i_cov_on = i_cov_on
          )
          fit_pois <- run_em_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            init_params = init_values_pois, 
            n_nodes = nodes,
            i_cov_on = i_cov_on
          )
        }
      } else if (length(i_cov_on) == 3) {
        init_values_pois <- get_start_values_poisson_with_cov(
          data = data,
          p_covariates = p_covariates,
          i_covariates = i_covariates,
          i_cov_on = c("alpha", "delta")
        )
        fit_pois <- run_em_poisson_with_cov(
          data = data,
          p_covariates = p_covariates,
          i_covariates = i_covariates,
          init_params = init_values_pois, 
          n_nodes = nodes,
          i_cov_on =  c("alpha", "delta")
        )
      } else if (length(i_cov_on) == 2) {
        if ("log_disp" %in% i_cov_on) {
          init_values_pois <- get_start_values_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            i_cov_on = i_cov_on[-which(i_cov_on == "log_disp")]
          )
          fit_pois <- run_em_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            init_params = init_values_pois, 
            n_nodes = nodes,
            i_cov_on =  i_cov_on[-which(i_cov_on == "log_disp")]
          )
        } else {
          init_values_pois <- get_start_values_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            i_cov_on = i_cov_on
          )
          fit_pois <- run_em_poisson_with_cov(
            data = data,
            p_covariates = p_covariates,
            i_covariates = i_covariates,
            init_params = init_values_pois, 
            n_nodes = nodes,
            i_cov_on =  i_cov_on
          )
        }
    }
    } else if (!is.null(p_covariates)) {
      init_values_pois <- get_start_values_poisson_with_cov(
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        i_cov_on = NULL
      )
      fit_pois <- run_em_poisson_with_cov(
        data = data,
        p_covariates = p_covariates,
        i_covariates = i_covariates,
        init_params = init_values_pois, 
        n_nodes = nodes,
        i_cov_on =  NULL
      )
    }
  }
  init_alphas <- fit_pois$params[grepl("alpha", names(fit_pois$params)) &
                                   !grepl("beta", names(fit_pois$params))]
  init_deltas <- fit_pois$params[grepl("delta", names(fit_pois$params)) &
                                   !grepl("beta", names(fit_pois$params))]
  
  if (!is.null(p_covariates)) {
    # we have a model with person covariates
    init_betas_p <- fit_pois$params[grepl("beta_p", names(fit_pois$params))]
    init_logdisps<-c()
    sim_abilities=rnorm(nsim)
    for (i in 1:ncol(data)) { 
      if (same_alpha) {
        mu <- exp(init_deltas[i] + init_alphas*sim_abilities) 
                    # + init_alphas*sum(t(init_betas_p * t(p_covariates))))
        # here alphas is a scalar because we have the constraint same alpha here
      } else {
        mu <- exp(init_deltas[i] + init_alphas[i]*sim_abilities) 
                  # +  init_alphas[i]*sum(t(init_betas_p * t(p_covariates))))
      }
      sim <- rpois(nsim, mu)
      init_logdisps[i] <- log((var(sim) / var(data[,i])))
    }
    
    start_values <- c(init_alphas, init_deltas, init_logdisps, init_betas_p)
    names(start_values) <- c(
      paste0("alpha", 1:length(init_alphas)),
      paste0("delta", 1:length(init_deltas)),
      paste0("log_disp", 1:length(init_logdisps)),
      paste0("beta_p", 1:length(init_betas_p))
    )
    
  } else if (!is.null(i_covariates)) {
    # we have a model with item covariates
    init_logdisps<-c()
    sim_abilities=rnorm(nsim)
    # distinguish between on which parameters we have item covariates
    if (length(i_cov_on) == 1) {
      if (i_cov_on == "delta") {
        init_betas_i <- fit_pois$params[grepl("beta_i", names(fit_pois$params))]
        for (i in 1:ncol(data)) {
          if (same_alpha) {
            mu <- exp(init_deltas + init_alphas*sim_abilities)
                        # + sum(t(init_betas_i * t(i_covariates))))
          } else {
            mu <- exp(init_deltas + init_alphas[i]*sim_abilities)
                        # + sum(t(init_betas_i * t(i_covariates))))
          }
          sim <- rpois(nsim, mu)
          init_logdisps[i] <- log((var(sim) / var(data[,i])))
        }
        start_values <- c(init_alphas, init_deltas, init_logdisps, init_betas_i)
        names(start_values) <- c(
          paste0("alpha", 1:length(init_alphas)),
          paste0("delta", 1:length(init_deltas)),
          paste0("log_disp", 1:length(init_logdisps)),
          paste0("beta_i", 1:length(init_betas_i))
        )
      } else if (i_cov_on == "alpha") {
        init_betas_i <- fit_pois$params[grepl("beta_i", names(fit_pois$params))]
        # we can't have the constraint of same alphas if we have covaraites on
        # alpha because the covariates have different values for the different
        # items implying different alphas
        for (i in 1:ncol(data)) {
          mu <- exp(init_deltas[i] + init_alphas*sim_abilities)
                     # + sim_abilities * sum(t(init_betas_i * t(i_covariates))))
          sim <- rpois(nsim, mu)
          init_logdisps[i] <- log((var(sim) / var(data[,i])))
        }
        start_values <- c(init_alphas, init_deltas, init_logdisps, init_betas_i)
        names(start_values) <- c(
          paste0("alpha", 1:length(init_alphas)),
          paste0("delta", 1:length(init_deltas)),
          paste0("log_disp", 1:length(init_logdisps)),
          paste0("beta_i", 1:length(init_betas_i))
        )
      } else if (i_cov_on == "log_disp") {
        for (i in 1:ncol(data)) {
          if (same_alpha) {
            mu <- exp(init_deltas[i] + init_alphas*sim_abilities)
          } else {
            mu <- exp(init_deltas[i] + init_alphas[i]*sim_abilities)
          }
          sim <- rpois(nsim, mu)
          init_logdisps[i] <- log((var(sim) / var(data[,i])))
        }
        init_logdisps <- mean(init_logdisps)
        # we need one log nu and then the covariate weights on nu
        predict_log_disp_df <- data.frame(
          count_var = log(apply(data, 2, var))
        )
        predict_log_disp_df <- as.data.frame(cbind(predict_log_disp_df, i_covariates))
        # i_covariates is a matrix with I columnds for I covaraites and M rows for the values
        # of the M items on those I covariates
        colnames(predict_log_disp_df)[-1] <- paste0("covar_", colnames(predict_log_disp_df)[-1])
        fit_log_disp <- lm(paste0("count_var ~", 
                                    paste(colnames(predict_log_disp_df)[-1], collapse = "+" )),
                             data = predict_log_disp_df)
        init_betas_i <- fit_log_disp$coefficients[-1]
        start_values <- c(init_alphas, init_deltas, init_logdisps, init_betas_i)
        names(start_values) <- c(
          paste0("alpha", 1:length(init_alphas)),
          paste0("delta", 1:length(init_deltas)),
          paste0("log_disp", 1:length(init_logdisps)),
          paste0("beta_i", 1:length(init_betas_i))
        )
      } 
    } else if (length(i_cov_on) == 3) {
      init_betas_i_delta <- fit_pois$params[grepl("beta_i_delta", names(fit_pois$params))]
      init_betas_i_alpha <- fit_pois$params[grepl("beta_i_alpha", names(fit_pois$params))]
      
      for (i in 1:ncol(data)) {
        # if we have covariates on all three item parameters, then we can't have any 
        # constraints but as a consequence of having covariates on all item parameters,
        # we have only scalars for init_alphas, init_deltas, and init_logdisps
        mu <- exp(init_deltas + init_alphas*sim_abilities)
        sim <- rpois(nsim, mu)
        init_logdisps[i] <- log((var(sim) / var(data[,i])))
      }
      
      init_logdisps <- mean(init_logdisps)
      # we need one log nu and then the covariate weights on nu
      predict_log_disp_df <- data.frame(
        count_var = log(apply(data, 2, var))
      )
      predict_log_disp_df <- as.data.frame(cbind(predict_log_disp_df, i_covariates))
      # i_covariates is a matrix with I columnds for I covaraites and M rows for the values
      # of the M items on those I covariates
      colnames(predict_log_disp_df)[-1] <- paste0("covar_", colnames(predict_log_disp_df)[-1])
      fit_log_disp <- lm(paste0("count_var ~", 
                                paste(colnames(predict_log_disp_df)[-1], collapse = "+" )),
                         data = predict_log_disp_df)
      init_betas_i_logdisp <- fit_log_disp$coefficients[-1]
      
      start_values <- c(init_alphas, init_deltas, init_logdisps, 
                        init_betas_i_alpha, init_betas_i_delta, init_betas_i_logdisp)
      names(start_values) <- c(
        paste0("alpha", 1:length(init_alphas)),
        paste0("delta", 1:length(init_deltas)),
        paste0("log_disp", 1:length(init_logdisps)),
        paste0("beta_i_alpha", 1:length(init_betas_i_alpha)),
        paste0("beta_i_delta", 1:length(init_betas_i_delta)),
        paste0("beta_i_log_disp", 1:length(init_betas_i_logdisp))
      )
    } else if (length(i_cov_on) == 2) {
      if ("log_disp" %in% i_cov_on) {
        # item covariates on alpha and log_disp or delta and log_disp
        if ("alpha" %in% i_cov_on) {
          # initial values for beta on alpha
          init_betas_i_alpha <- fit_pois$params[grepl("beta_i", names(fit_pois$params))]
          # the b weight won't be named after alpha here because in this poisson model, we only
          # have covariates on alpha then here
          # initial values for intercept on log_disp
          for (i in 1:ncol(data)) {
            mu <- exp(init_deltas[i] + init_alphas*sim_abilities)
            # + sim_abilities * sum(t(init_betas_i * t(i_covariates))))
            sim <- rpois(nsim, mu)
            init_logdisps[i] <- log((var(sim) / var(data[,i])))
          }
          init_logdisps <- mean(init_logdisps)
          # initial values for beta on log_disp
          predict_log_disp_df <- data.frame(
            count_var = log(apply(data, 2, var))
          )
          predict_log_disp_df <- as.data.frame(cbind(predict_log_disp_df, i_covariates))
          # i_covariates is a matrix with I columnds for I covaraites and M rows for the values
          # of the M items on those I covariates
          colnames(predict_log_disp_df)[-1] <- paste0("covar_", colnames(predict_log_disp_df)[-1])
          fit_log_disp <- lm(paste0("count_var ~", 
                                    paste(colnames(predict_log_disp_df)[-1], collapse = "+" )),
                             data = predict_log_disp_df)
          init_betas_i_logdisp <- fit_log_disp$coefficients[-1]
          # prepare start values for ouput
          start_values <- c(init_alphas, init_deltas, init_logdisps, 
                            init_betas_i_alpha, init_betas_i_logdisp)
          names(start_values) <- c(
            paste0("alpha", 1:length(init_alphas)),
            paste0("delta", 1:length(init_deltas)),
            paste0("log_disp", 1:length(init_logdisps)),
            paste0("beta_i_alpha", 1:length(init_betas_i_alpha)),
            paste0("beta_i_log_disp", 1:length(init_betas_i_logdisp))
          )
        } else {
          # then we have the pairing delta and log_disp
          # initial values for bet on delta
          init_betas_i_delta <- fit_pois$params[grepl("beta_i", names(fit_pois$params))]
          # the b weight won't be named after delta here because in this poisson model, we only
          # have covariates on delta then here
          # initial values for intercept on log disp
          for (i in 1:ncol(data)) {
            if (same_alpha) {
              mu <- exp(init_deltas + init_alphas*sim_abilities)
              # + sim_abilities * sum(t(init_betas_i * t(i_covariates))))
            } else {
              mu <- exp(init_deltas + init_alphas[i]*sim_abilities)
              # + sim_abilities * sum(t(init_betas_i * t(i_covariates))))
            }
            sim <- rpois(nsim, mu)
            init_logdisps[i] <- log((var(sim) / var(data[,i])))
          }
          init_logdisps <- mean(init_logdisps)
          # initial values for beta on log_disp
          predict_log_disp_df <- data.frame(
            count_var = log(apply(data, 2, var))
          )
          predict_log_disp_df <- as.data.frame(cbind(predict_log_disp_df, i_covariates))
          # i_covariates is a matrix with I columnds for I covaraites and M rows for the values
          # of the M items on those I covariates
          colnames(predict_log_disp_df)[-1] <- paste0("covar_", colnames(predict_log_disp_df)[-1])
          fit_log_disp <- lm(paste0("count_var ~", 
                                    paste(colnames(predict_log_disp_df)[-1], collapse = "+" )),
                             data = predict_log_disp_df)
          init_betas_i_logdisp <- fit_log_disp$coefficients[-1]
          # prepare start values for ouput
          start_values <- c(init_alphas, init_deltas, init_logdisps, 
                            init_betas_i_delta, init_betas_i_logdisp)
          names(start_values) <- c(
            paste0("alpha", 1:length(init_alphas)),
            paste0("delta", 1:length(init_deltas)),
            paste0("log_disp", 1:length(init_logdisps)),
            paste0("beta_i_delta", 1:length(init_betas_i_delta)),
            paste0("beta_i_log_disp", 1:length(init_betas_i_logdisp))
          )
        }
      } else {
        # item covariates on alpha and delta
        init_betas_i_delta <- fit_pois$params[grepl("beta_i_delta", names(fit_pois$params))]
        init_betas_i_alpha <- fit_pois$params[grepl("beta_i_alpha", names(fit_pois$params))]
        # start values for the log_disps
        for (i in 1:ncol(data)) {
          mu <- exp(init_deltas + init_alphas*sim_abilities)
          # + sim_abilities * sum(t(init_betas_i * t(i_covariates))))
          sim <- rpois(nsim, mu)
          init_logdisps[i] <- log((var(sim) / var(data[,i])))
        }
        # prepare start values for ouput
        start_values <- c(init_alphas, init_deltas, init_logdisps, 
                          init_betas_i_alpha, init_betas_i_delta)
        names(start_values) <- c(
          paste0("alpha", 1:length(init_alphas)),
          paste0("delta", 1:length(init_deltas)),
          paste0("log_disp", 1:length(init_logdisps)),
          paste0("beta_i_alpha", 1:length(init_betas_i_alpha)),
          paste0("beta_i_delta", 1:length(init_betas_i_delta))
        )
      }
    }
  }
  
  return(start_values)
}

