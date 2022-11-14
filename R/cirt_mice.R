#' Model fitting function for count data IRT models using \link[countirt]{cirt} on mids objects (as returned after performing multiple imputations with \link[mice]{mice}).
#' 
#' @param model A string specifying the model. See \link[countirt]{cirt} for more information.
#' @param mids_object A mids object as obtained after performing multiple imputations on a dataset using \link[mice]{mice}.
#' @param family A string indicating the count data family, can be either "cmp" or "poisson".
#' @param item_offset Either a scalar (for same offset for all items) or a vector of the same length as the number of items with an offset to be added to the prediction term (on log scale).
#' @param person_offset Either a scalar (for same offset for all persons) or a vector of the same length as the number of persons with an offset to be added to the prediction term (on log scale).
#' @param data_long A boolean. Indicates whether data is in long format. If FALSE, expects data in wide format. Defaults to FALSE.
#' @param person_id A character string. Name of the column with person id in long format data frame. Only necessary if data_long = TRUE.
#' @param control A list providing control parameters for \link[countirt]{cirt}.
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
#' @importFrom mice complete
#' @useDynLib countirt, .registration=TRUE
#' @export
cirt_mice <- function(model, mids_object, family,
                 item_offset = NULL,
                 person_offset = NULL,
                 data_long = FALSE,
                 person_id = NULL,
                 control = list(
                   n_nodes = 121,
                   thres = Inf, prob = 0, init_disp_one = TRUE,
                   maxiter = 1000, 
                   convtol = 1e-5, ctol_maxstep = 1e-8,
                   m_method = "nleqslv", convcrit = "marglik"
                 )) {
  
  # fit the 2PCMPM to each imputed data set
  imputed_fits <- vector(mode = "list", length = mids_object$m)
  for (i in 1:mids_object$m) {
    cat(sprintf('\rImputed data set: %d', i))
    df <- complete(mids_object, i)
    imputed_fits[[i]] <- cirt(
      model = model,
      data = df,
      family = family,
      item_offset = item_offset,
      person_offset = person_offset,
      data_long = data_long,
      person_id = person_id,
      control = control
    )
    imputed_fits[[i]] <- add_inference(imputed_fits[[i]])
  }
  
  # pool effect estimates and standard errors
  # (under assumption that parameter estimates are normally distributed)
  params_df <- data.frame(
    imp = 1:mids_object$m
  )
  for (i in 1:length(imputed_fits[[1]]$fit$params)) {
    params_df[[names(imputed_fits[[1]]$fit$params)[i]]] <- NA
    params_df[[paste0("se_", names(imputed_fits[[1]]$fit$params)[i])]] <- NA
  }
  for (i in 1:mids_object$m) {
    params_df[i,names(imputed_fits[[1]]$fit$params)] <- imputed_fits[[i]]$fit$params
    params_df[i,paste0("se_", names(imputed_fits[[1]]$fit$params))] <- imputed_fits[[i]]$fit_ses$se
  }
  # compute average effect estimates across imputed data sets
  pooled_params <- apply(params_df[,names(imputed_fits[[1]]$fit$params)], 2, mean)
  names(pooled_params) <- names(imputed_fits[[1]]$fit$params)
  # compute average standard errors across imputed data sets
  # within imputation variance
  avg_within_var <- apply(params_df[,paste0("se_", names(imputed_fits[[1]]$fit$params))], 
                          2, 
                          function(x){mean(x^2)})
  names(avg_within_var) <- names(imputed_fits[[1]]$fit$params)
  # between imputation variance
  avg_btw_var <- apply(params_df[,names(imputed_fits[[1]]$fit$params)], 2, var)
  names(avg_btw_var) <- names(imputed_fits[[1]]$fit$params)
  # total variance = v_within + v_between + v_between / m
  pooled_var <- avg_within_var + avg_btw_var  + (avg_btw_var/mids_object$m)
  # pooled ses
  pooled_se <- sqrt(pooled_var)
  
  # compute pooled wald tests and p-values based on t distribution
  pooled_wald <- pooled_params / pooled_se
  # compute df for pooled_wald sampling t distribution
  # quantile: 1-alpha/2
  riv <- (avg_btw_var + (avg_btw_var/mids_object$m)) / avg_within_var
  # riv is the relative increase in variance due to nonrespnse (RIV)
  df_old <- (mids_object$m-1) * (1+(1/riv))^2
  lambda <- (avg_btw_var + (avg_btw_var/mids_object$m)) / pooled_var
  n <- nrow(complete(mids_object, 1))
  k <- length(imputed_fits[[1]]$fit$params)
  df_obs <- (((n - k) + 1) / ((n - k) + 3)) * (n - k)*(1 - lambda)
  df_adj <- (df_old * df_obs) / (df_old + df_obs)
  pvalues <- 2*(1 - pt(abs(pooled_wald), df_adj))
  
  # confidence intervals for pooled parameter estimates
  prob <- 0.95
  CI_quantiles <- -qt((1-prob)/2, df_adj)
  pooled_CI_lower <- pooled_params - CI_quantiles*pooled_se
  pooled_CI_upper <- pooled_params + CI_quantiles*pooled_se
  
  # prepare output
  out <- list(
    pooled_params = pooled_params,
    pooled_se = pooled_se,
    pooled_CI_lower = pooled_CI_lower,
    pooled_CI_upper = pooled_CI_upper,
    pooled_wald = pooled_wald,
    df_adj = df_adj,
    pvalues = pvalues,
    imputed_fits = imputed_fits
  )
  
  return(out)
}  
  
  
  
  
  
  