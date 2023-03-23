#' Model fitting function for count data IRT models.
#' 
#' @param model A string specifying the model. 
#' 
#' - For a data frame in wide format (default; for the 2PCMP and the CLRM): You specify the items that belong to the factor by writing "theta =~ item1 + item2 ...;" where you replace "itemj" with correct name in your data frame in wide format (note that currently, only one-dimensional models are possible but in the future, I might also implement multi-dimensional models). You can specify constraints (as long as they are implemented, see below) on the item parameters. For constraints on the slope parameters, you can currently contrain them to be equal across item ("alphas ~ 1;") or fix them to specific values ("alpha1 =: 0.5 alpha2 =: 0.3 ...;" But please note that currently, either all slopes must be fixed to specific values or all are estimated. This will be extended in the future to allow for some items to have fixed slopes and others to have freely estimated slopes.). The same is possible for the intercepts delta (with "deltas ~ 1;"; but note that fixing intercepts to certain values has not yet been implemented) and the (log) dispersions log_nu (with "log_nus ~ 1;" and "log_nu1 =: 0.2 log_nu =: 1.1 ...;"). If you wish to specify a CLRM, you can additionally include covariates on the latent ability with "thetas ~ 1 + C1 + C2 + ...;". The covariates must be named as in the data frame. The intercept will be a random intercept (with latent variance fixed to 1 and a mean of 0). Please note that at the moment, I would only recommend  using categorical covariates as computation for continuous covariates is very slow (albeit implemented). You can specify the formula in the same way for categorical covariates, just make sure that all categorical covariates in your data frame are factors.
#' 
#' - For a data frame in long format (note that you must specify the argument `long_data = TRUE`; for the DRTM): Specify the factor formula as "theta =~ counts(itemid::item1) + counts(itemid::item2) ...;" where you replace "count" with the name of the column with responses in your data frame in long format and "itemj" with the correct item identifier as in your data frame in long format and "itemid" with the name of the column in your data frame in long format that contains the item identifier. You need to use the :: between identifier column and respective item identifier. You can specify item parameter constraints as described above. To include covariates on any of the three item parameters, specify "alphas ~ 1 + C1 + C2 ...;" (or "deltas ~ ...;" or "log_nus ~ ...;" respectively), where you replace Ci with the names of the covariates in your data frame in long format. You can include covariates on any subset of parameters and constraints (as described above) on the remaining item parameters. Note that it is not possible to include item parameters on constrained parameters though. At the moment, you can choose which item parameters you want to include the covariates on but if you want to include covariates on more than one item parameter, they have to be the same across parameters at the moment. In the future, I will probably allow more flexibility in this regard.
#' 
#' Please note that the model specification expects you to use the name for the parameters as explained here. Your items must have the same name as in the data frame. If you just specify the factor formula, the full 2PCMPM with no covariates will be estimated. Each line of specification must start with the correct parameter name as explained above and must end with ; . Note that for the model specification to be successfully passed, the variable names in your data frame must not contain any '(' nor any '::'.
#' @param data A data frame either in long or in wide format. For the 2PCMP and the CLRM, wide format is required (which is the default expected format). For the DRTM, long format is required (for which you have to specify `long_format = TRUE`).
#' @param family A string indicating the count data family, can be either "cmp" or "poisson".
#' @param item_offset Either a scalar (for same offset for all items) or a vector of the same length as the number of items with an offset to be added to the prediction term (on log scale).
#' @param person_offset Either a scalar (for same offset for all persons) or a vector of the same length as the number of persons with an offset to be added to the prediction term (on log scale).
#' @param data_long A boolean. Indicates whether data is in long format. If FALSE, expects data in wide format. Defaults to FALSE.
#' @param person_id A character string. Name of the column with person id in long format data frame. Only necessary if data_long = TRUE.
#' @param stand_errors A boolean. Indicates whether standard errors for model parameters should be estimated. Defaults to FALSE.
#' @param control A list providing control parameters for the estimation.
#' 
#' @import Rcpp
#' @import RcppGSL
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom fastGHQuad ghQuad
#' @importFrom nleqslv nleqslv
#' @importFrom rootSolve gradient
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr enquo
#' @importFrom turboEM turboem
#' @useDynLib countirt, .registration=TRUE
#' @export
cirt <- function(model, data, family,
                 item_offset = NULL,
                 person_offset = NULL,
                 data_long = FALSE,
                 person_id = NULL,
                 stand_errors = FALSE,
                 accelerate = "none",
                 control = list(
                   n_nodes = 121,
                   thres = Inf, prob = 0, init_disp_one = TRUE,
                   maxiter = 1000, 
                   convtol = 1e-5, ctol_maxstep = 1e-8,
                   m_method = "nleqslv", convcrit = "marglik"
                 )) {
  # TODO checks and data prep
  # TODO implement some proper error catching and meaningful error messages
  # TODO add verbose argument and only print iterations if TRUE
  # TODO print nicer when printing iterations
  
  # TODO work on accelerator implementation (had to take it out for now)
  # add argument back in: accelerate_method = replicate(length(accelerate),list())
  # #' #' @param accelerate A string. Can bei either "none" for no acceleration or "squarem"
  # , "pem", "qn" for accelerators as implemented in turboEM. Control parameters can be set 
  # via the accelerate_method element in the control list. Compare turboem documentation for
  #how to specify control parameters for the respective methods (accelerate_method will be 
  #passed to control.method element in turbem).
  
  data <- as.data.frame(data)
  
  # extract model info from model specification -------------------------
  model_list <- parse_model(
    model = model, data = data, 
    data_long = data_long, person_id = person_id,
    family = family
    )
  
  # model list checks ---------------------------------------------------
  
  if (!is.null(item_offset)) {
    # if we have item_offsets, check that they make sense
    if (length(item_offset) == 1) {
      # if we just one scalar provided as item intercept, then we need to turn
      # that into a vector of length n_items for all subsequent functions
      item_offset <- rep(item_offset, ncol(model_list$item_data))
    } else {
      # if not just one scalar, make sure that we have as many entries in
      # item_offset as we have items, otherwise throw an error
      if (length(item_offset) != ncol(model_list$item_data)) {
        stop("Argument item_offset must either be a scalar to be used as the offset for all items or must be the same length as the number of items used in the model. Currently, neither is the case. Please check and adjust the item_offset argument.")
      }
    }
  }
  
  if (!is.null(person_offset)) {
    # if we have personoffsets, check that they make sense
    if (length(person_offset) == 1) {
      # if we just one scalar provided as person offser, then we need to turn
      # that into a vector of length n_persons for all subsequent functions
      person_offset <- rep(person_offset, nrow(model_list$item_data))
    } else {
      # if not just one scalar, make sure that we have as many entries in
      # person_offset as we have persons, otherwise throw an error
      if (length(person_offset) != nrow(model_list$item_data)) {
        stop("Argument person_offset must either be a scalar to be used as the offset for all persons or must be the same length as the number of persons in the sample. Currently, neither is the case. Please check and adjust the person_offset argument.")
      }
    }
  }
  
  if (!is.null(model_list$p_covariates) & !is.null(model_list$i_covariates)) {
    stop("Currently, countirt does not support person and item covariates at the same time. This might be implemented in the future, but for now, please specify a model with either only item or only person covariates.")
  }
  
  if (!is.null(model_list$i_covariates)) {
    if (family == "poisson" & ("log_disp" %in% model_list$i_cov_on)) {
      # check that no item covariates on the dispersion were specified when a Poisson model is fitted 
      # as we don't have dispersions in the poisson model
      stop("You specified item covariates on the dispersion parameters, but specified family = \"poisson\". There are no dispersion parameters in a Poisson model. Please either remove the item covariates on the dispersion, or change the family to \"cmp\".")
    }
  }
    
  

  # TODO checks that we don't specify constraints on parameters with covariate
  
  # TODO wenn ich hier das i_cov_on argument einbaue, so dass ich item kovariaten mit
  # einbauen kann, dann sollte ich checken, dass wenn fix_alphas = TRUE, nur
  # i_cov_on = "delta" ist; wenn nicht, das so setzen und eine warnung ausgeben, dass
  # man bei fixierten alphas dann nur item kovaraiten auf den deltas haben kann
  # (im poisson fall); analog auch im cmp fall unter beruecksichtigung von nu
  # (same with same_alpha; and analgously fix_disps and same_disp)
  
  # TODO check for i_cov_on %in% c("alpha", "delta", "log_nu")
  # also dass wir keinen unsinn haben als angegebenes element
  # und auch entspr. fuer poisson
  
  # handling missings ---------------------------------------------------
  
  if (any(is.na(model_list$item_data))) {
    if (!is.null(model_list$p_covariates)) {
      if (any(is.na(model_list$p_covariates))) {
        # exclude rows with missings on either items or on person covariates
        # so that we have full obs for all persons included in a model
        # with person covariates
        missings_item_df <- which(is.na(model_list$item_data))
        missings_pcov_df <- which(is.na(model_list$p_covariates))
        missings <- unique(c(missings_item_df,missings_pcov_df))
        model_list$item_data <- model_list$item_data[-missings,]
        model_list$p_covariates <- model_list$p_covariates[-missings,]
        warning("Rows with missing values on either items or person covariates or both have been removed.")
      } else {
        # no missings on person covariates, just on data
        model_list$item_data <- na.omit(model_list$item_data)
        warning("Rows with missing values on item data have been removed.")
      }
    } else { # we don't have person covariates
      model_list$item_data <- na.omit(model_list$item_data)
      warning("Rows with missing values have been removed.")
    }
  }
  
  # check if we have missings just on person covariates even if we have
  # complete obs only on item data
  if (!is.null(model_list$p_covariates)) {
    if (any(is.na(model_list$p_covariates))) {
      # exclude rows with missings on either items or on person covariates
      # so that we have full obs for all persons included in a model
      # with person covariates
      missings <- which(is.na(model_list$p_covariates))
      model_list$item_data <- model_list$item_data[-missings,]
      model_list$p_covariates <- model_list$p_covariates[-missings,]
      warning("Rows with missing values on person covariates been removed.")
    }
  }
  
  # check for missings on item covariates and break if we have those as we
  # the drtm can't handle that, or otherwise we would have to remove the entire
  # item and i wouldn't wanna make that choice for the user.
  if (!is.null(model_list$i_covariates)) {
    if (any(is.na(model_list$i_covariates))) {
      stop("There are missing values on at least one of the item covariates. Each item covariate must have a value for each item. Otherwise please remove the item covariate(s) with the missings or remove the item for which covariates have missings.")
    }
  }
  
  # model fitting -----------------------------------------------------------------
  
  if (is.null(model_list$p_covariates) & is.null(model_list$i_covariates)) {
    # 2pcmp model
    if (family == "cmp") {
      print("Start determining start values.")
      
      if (!is.null(model_list$fixed_log_disps)) {
        fixed_disps <- exp(model_list$fixed_log_disps)
      } else {
        fixed_disps <- NULL
      }
      
      # TODO init_disp_one argument entfernen; auch aus der get_start_value funktion
      start_values <- get_start_values(
        data = model_list$item_data,
        init_disp_one = control$init_disp_one,
        fix_disps = fixed_disps, 
        fix_alphas = model_list$fixed_alphas,
        same_disps = model_list$equal_log_disps, 
        same_alpha = model_list$equal_alphas,
        item_offset = item_offset,
        person_offset = person_offset
      )
      
      print("Start model fitting. This will take a little bit of time.")
      #if (accelerate == "none") {
        fit <- run_newem(
          init_params = start_values,
          data = model_list$item_data, 
          n_nodes = control$n_nodes, 
          fix_disps = fixed_disps,
          fix_alphas = model_list$fixed_alphas,
          same_disps = model_list$equal_log_disps, 
          same_alphas = model_list$equal_alphas,
          item_offset = item_offset,
          person_offset = person_offset,
          thres = control$thres,
          prob = control$prob,
          maxiter = control$maxiter, 
          convtol = control$convtol, 
          ctol_maxstep = control$ctol_maxstep,
          m_method = control$m_method, 
          convcrit = control$convcrit
        )
      # } else {
      #   weights_and_nodes <- quad_rule(control$n_nodes, 
      #                                 thres = control$thres, 
      #                                 prob = control$prob)
      #   if (is.null(item_offset)) {
      #     item_offset <- rep(0, ncol(data))
      #   }
      #   if (is.null(person_offset)) {
      #     person_offset <- rep(0, nrow(data))
      #   }
      #   fit <- turboem(
      #     par = start_values,
      #     fixtptfn = newem_em_cycle,
      #     objfn = neg_marg_ll,
      #     method = accelerate,
      #     control.method = control$accelerate_method,
      #     # arguments for fixtptfn and objfn
      #     data = model_list$item_data,
      #     weights_and_nodes = weights_and_nodes,
      #     ctol_maxstep = control$ctol_maxstep, 
      #     m_method = control$m_method,
      #     fix_disps = fixed_disps, 
      #     fix_alphas = model_list$fixed_alphas,
      #     same_disps = model_list$equal_log_disps, 
      #     same_alphas = model_list$equal_alphas, 
      #     item_offset = item_offset,
      #     person_offset = person_offset,
      #     family = "cmp"
      #   )
        # TODO 
        # this is fine for outputting out of cirt, i just need to support this in summary
        # i have to think about whether i am shaping this here (which is probably nicer)
        # or in summary
      #}

    } else if (family == "poisson") {
      start_values <- get_start_values_pois(
        data = model_list$item_data,
        same_alpha = model_list$equal_alphas,
        fix_alphas = model_list$fixed_alphas,
        item_offset = item_offset,
        person_offset = person_offset
      )
      
      print("Start model fitting.")
      #if (accelerate == "none") {
        fit <- run_em_poisson(
          init_params = start_values, 
          data = model_list$item_data,
          n_nodes = control$n_nodes, 
          fix_alphas = model_list$fixed_alphas,
          same_alpha = model_list$equal_alphas,
          item_offset = item_offset,
          person_offset = person_offset,
          thres = control$thres,
          prob = control$prob,
          maxiter = control$maxiter, 
          convtol = control$convtol, 
          ctol_maxstep = control$ctol_maxstep,
          convcrit = control$convcrit
        )
      # } else {
      #   weights_and_nodes<- quad_rule(control$n_nodes, 
      #                                 thres = control$thres, 
      #                                 prob = control$prob)
      #   if (is.null(item_offset)) {
      #     item_offset <- rep(0, ncol(data))
      #   }
      #   if (is.null(person_offset)) {
      #     person_offset <- rep(0, nrow(data))
      #   }
      #   fit <- turboem(
      #     par = start_values,
      #     fixtptfn = em_cycle_poisson,
      #     objfn = neg_marg_ll,
      #     method = accelerate,
      #     control.method = control$accelerate_method,
      #     # arguments for fixtptfn and objfn
      #     data = model_list$item_data,
      #     weights_and_nodes = weights_and_nodes,
      #     ctol_maxstep = control$ctol_maxstep, 
      #     m_method = mcontrol$m_method,
      #     fix_alphas = model_list$fixed_alphas,
      #     same_alphas = model_list$equal_alphas, 
      #     item_offset = item_offset,
      #     person_offset = person_offset,
      #     family = "poisson"
      #   )
        # TODO 
        # this is fine for outputting out of cirt, i just need to support this in summary
        # i have to think about whether i am shaping this here (which is probably nicer)
        # or in summary
      #}
      
    } else {
      stop("Invalid family specified. Please see documentations for available families.")
    }
  } else {
    # 2pcmp model with covariates
    # TODO add offset option here for covariate models
    if (family == "cmp") {
      print("Start determining start values.")
      # TODO think about whether I am dealing correctly here with all possible constraints,
      # what about fixed alphas and/or disps
      # TODO i think i still need to implement fix alphas into start values for covariate models
      start_values <- get_start_values_cmp_with_cov(
        data = model_list$item_data,
        p_covariates = model_list$p_covariates,
        i_covariates = model_list$i_covariates,
        nodes = control$n_nodes,
        same_alpha = model_list$equal_alphas,
        i_cov_on = model_list$i_cov_on,
        which_i_cov = model_list$which_i_cov, 
        item_offset = item_offset
      )
      
      if (!is.null(model_list$fixed_log_disps)) {
        fixed_disps <- exp(model_list$fixed_log_disps)
      } else {
        fixed_disps <- NULL
      }
      
      print("Start model fitting. This will take a little bit of time.")
      fit <- run_em_cmp_with_cov(
        data = model_list$item_data, 
        init_params = start_values, 
        n_nodes = control$n_nodes, 
        p_covariates = model_list$p_covariates,
        i_covariates = model_list$i_covariates,
        i_cov_on = model_list$i_cov_on,
        which_i_cov = model_list$which_i_cov, 
        p_cov_cat = model_list$p_cov_cat,
        num_levels_p_cov = model_list$p_cov_levels,
        fix_disps = fixed_disps,
        fix_alphas = model_list$fixed_alphas,
        same_disps = model_list$equal_log_disps, 
        same_alphas = model_list$equal_alphas,
        item_offset = item_offset,
        thres = control$thres,
        prob = control$prob,
        maxiter = control$maxiter, 
        convtol = control$convtol, 
        ctol_maxstep = control$ctol_maxstep,
        m_method = control$m_method, 
        convcrit = control$convcrit
      )
      
    } else if (family == "poisson") {
      
      start_values <- get_start_values_poisson_with_cov(
        data = model_list$item_data,
        p_covariates = model_list$p_covariates,
        i_covariates = model_list$i_covariates,
        same_alpha = model_list$equal_alphas,
        i_cov_on = model_list$i_cov_on,
        which_i_cov = model_list$which_i_cov,
        item_offset = item_offset
      )
      
      print("Start model fitting.")
      fit <- run_em_poisson_with_cov(
        data = model_list$item_data, 
        init_params = start_values, 
        n_nodes = control$n_nodes, 
        p_covariates = model_list$p_covariates,
        i_covariates = model_list$i_covariates,
        i_cov_on = model_list$i_cov_on,
        which_i_cov = model_list$which_i_cov, 
        fix_alphas = model_list$fixed_alphas,
        same_alpha = model_list$equal_alphas,
        item_offset = item_offset,
        thres = control$thres,
        prob = control$prob,
        maxiter = control$maxiter, 
        convtol = control$convtol, 
        ctol_maxstep = control$ctol_maxstep,
        convcrit = control$convcrit
      )
      
    } else {
      stop("Invalid family specified. Please see documentations for available families.")
    }
  }
  
  # prepare object for returning
  out <- list(
    family = family,
    model = model_list,
    fit = fit,
    fit_ses = NULL, # we can add ses with add_inference
    start_values = start_values,
    accelerate = accelerate,
    control = control
  )
  class(out) <- "cirtfit"
  return(out)
}

# TODO check that all possible constraints are implemented and are also
# implemented all throughout, i.e., so that we can have standard errors for all
# possible models we can fit, including which_i_cov and poisson variants etc.

#' Function that computes standard errors, z- and p-values for Wald tests, and CI.
#' 
#' @param model A cirt model. As returned by the cirt function.
#' @param prob Probability for the CI. Defaults to 0.95
#' 
#' @import Rcpp
#' @import RcppGSL
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom fastGHQuad ghQuad
#' @importFrom nleqslv nleqslv
#' @importFrom rootSolve gradient
#' @useDynLib countirt, .registration=TRUE
#' @export
add_inference <- function(model, prob = 0.95) {
  
  if (!is.null(model$model$fixed_log_disps)) {
    fixed_disps <- exp(model$model$fixed_log_disps)
  } else {
    fixed_disps <- NULL
  }
  
  if (is.null(model$model$p_covariates) & is.null(model$model$i_covariates)) {
    # 2pcmpm without covariates
    if (model$family == "cmp") {
      vcov <- compute_vcov(
        item_params = model$fit$params,
        weights_and_nodes = quad_rule(model$control$n_nodes),
        data = model$model$item_data,
        same_alphas = model$model$equal_alphas, 
        same_disps = model$model$equal_log_disps,
        fix_alphas = model$model$fixed_alphas, 
        fix_disps = model$model$fixed_log_disps,
        item_offset = model$fit$item_offset
      )
    } else if (model$family == "poisson") {
      vcov <- compute_vcov_poisson(
        item_params = model$fit$params,
        weights_and_nodes = quad_rule(model$control$n_nodes),
        data = model$model$item_data,
        same_alphas = model$model$equal_alphas, 
        fix_alphas = model$model$fixed_alphas, 
        item_offset = model$fit$item_offset
      ) 
    }
  } else {
    # with covariates, so drtm or clrm
    if (model$family == "cmp") {
      if (isTRUE(model$model$p_cov_cat)) {
        resp_patterns_matrix <- make_resp_patterns_mat(
          lapply(model$model$p_cov_levels, get_resp_patterns_pcov_cat), 
          prod(model$model$p_cov_levels), 
          model$model$p_cov_levels,
          item_offset = model$fit$item_offset
        )
      } else {
        resp_patterns_matrix <- NULL
      }
      
      vcov <- compute_vcov_with_cov(
        item_params = model$fit$params,
        weights_and_nodes = quad_rule(model$control$n_nodes),
        data = model$model$item_data,
        p_covariates = model$model$p_covariates,
        i_covariates = model$model$i_covariates,
        i_cov_on = model$model$i_cov_on,
        which_i_cov = model$model$which_i_cov,
        p_cov_cat = model$model$p_cov_cat,
        resp_patterns_matrix = resp_patterns_matrix,
        same_alphas = model$model$equal_alphas, 
        same_disps = model$model$equal_log_disps,
        fix_alphas = model$model$fixed_alphas, 
        fix_disps = fixed_disps,
        item_offset = model$fit$item_offset
      )
    } else if (model$family == "poisson") {
      vcov <- compute_vcov_poisson_with_cov(
        item_params = model$fit$params,
        weights_and_nodes = quad_rule(model$control$n_nodes),
        data = model$model$item_data,
        p_covariates = model$model$p_covariates,
        i_covariates = model$model$i_covariates,
        i_cov_on = model$model$i_cov_on,
        same_alphas = model$model$equal_alphas, 
        fix_alphas = model$model$fixed_alphas,
        item_offset = model$fit$item_offset
      )
    }
  }
  
  CI_quantile <- -qnorm((1-prob)/2)
  se <- se_from_vcov(vcov)
  CI_lower <- model$fit$params - CI_quantile*se
  CI_upper <- model$fit$params + CI_quantile*se
  z_value <- model$fit$params/se
  p_value <- 2*(1-pnorm(abs(z_value)))
  
  inf_list <- list(
    se = se,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    z_value = z_value,
    p_value = p_value
  )
  
  out <- model
  out$fit_ses <- inf_list
  class(out) <- "cirtfit"
  return(out)
}

# TODO generic summary schreiben, wenn ich die auf
# cirt objekt aufrufe ohne se teil (is.null), weil ich den in cirt immer auf null setzen,
# warnung schreiben mit hinweis, dass ich das hinzufuegen kann mit add_inference()


# TODO wrapper for ability estimation



