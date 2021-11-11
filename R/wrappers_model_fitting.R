#' Model fitting function for count data IRT models
#' 
#' @param data a data matrix, each column must correspond to one item, each row to an observational unit (usually a person)
#' @param family a string indicating the count data family, can be either "cmp" or "poisson"
#' @param n_nodes an integer, number of quadrature nodes, defaults to 121, for the 2CMP model, no less than 100 quadrature nodes are recommended
#' @param stand_errors a boolean, indicates whether standard errors for model parameters should be estimated, defaults to FALSE
#' @param constraints a list, indicating the constraints for the model, note that at the current moment, the possible constraints are mutually exclusive (this functionality will be extended in the future). 
#' @param control a list, providing control parameters for the estimation
#' 
#' @import Rcpp
#' @import RcppGSL
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom fastGHQuad ghQuad
#' @importFrom nleqslv nleqslv
#' @importFrom rootSolve gradient
#' @useDynLib countirt, .registration=TRUE
#' @export
cirt <- function(data, family, n_nodes = 121, stand_errors = FALSE,
                 constraints = list(
                   fix_disps = NULL, fix_alphas = NULL,
                   same_disps = FALSE, same_alphas = FALSE
                 ),
                 control = list(
                   thres = Inf, prob = 0, init_disp_one = TRUE,
                   maxiter = 1000, 
                   convtol = 1e-5, ctol_maxstep = 1e-8,
                   m_method = "nleqslv", convcrit = "marglik"
                 )) {
  # TODO checks and data prep
  # TODO implement some proper error catching and meaningful error messages
  
  # TODO add some formula syntax to cirt
  # and think about what makes sense for being able to add item and person covariates
  
  if (any(is.na(data))) {
    # TODO remove rows with NAs and print warning that they were removed
  }
  
  # TODO wenn ich hier das i_cov_on argument einbaue, so dass ich item kovariaten mit
  # einbauen kann, dann sollte ich checken, dass wenn fix_alphas = TRUE, nur
  # i_cov_on = "delta" ist; wenn nicht, das so setzen und eine warnung ausgeben, dass
  # man bei fixierten alphas dann nur item kovaraiten auf den deltas haben kann
  # (im poisson fall); analog auch im cmp fall unter beruecksichtigung von nu
  # (same with same_alpha; and analgously fix_disps and same_disp)
  
  # TODO incorporate check that we only have item or person parameters as
  # we can only do one or the other atm
  
  # TODO check for i_cov_on %in% c("alpha", "delta", "log_nu")
  # also dass wir keinen unsinn haben als angegebenes element
  # und auch entspr. fuer poisson
  
  if (family == "cmp") {
    print("Start determining start values.")
    start_values <- get_start_values(
      data = data, init_disp_one = control$init_disp_one
    )
    
    print("Start model fitting. This will at least take a little bit of time.")
    fit <- run_newem(
      data = data, 
      init_params = start_values, 
      n_nodes = n_nodes, 
      fix_disps = constraints$fix_disps,
      fix_alphas = constraints$fix_alphas,
      same_disps = constraints$same_disps, 
      same_alphas = constraints$same_alphas,
      thres = control$thres,
      prob = control$prob,
      maxiter = control$maxiter, 
      convtol = control$convtol, 
      ctol_maxstep = control$ctol_maxstep,
      m_method = control$m_method, 
      convcrit = control$convcrit
    )
  
    if (stand_errors) {
      fit_vcov <- compute_vcov(fitparams, quad_rule(n_nodes), data)
      fit_ses <- se_from_vcov(fit_vcov)
    } else {
      fit_vcov <- NA
      fit_ses <- NA
    }
  } else if (family == "poisson") {
    start_values <- get_start_values_pois(
      data = data
    )
    
    print("Start model fitting.")
    fit <- run_em_poisson(
      data = data, 
      init_params = start_values, 
      n_nodes = n_nodes, 
      fix_alphas = constraints$fix_alphas,
      same_alpha = constraints$same_alphas,
      thres = control$thres,
      prob = control$prob,
      maxiter = control$maxiter, 
      convtol = control$convtol, 
      ctol_maxstep = control$ctol_maxstep,
      convcrit = control$convcrit
    )
    
    if (stand_errors) {
      fit_vcov <- compute_vcov_poisson(fitparams, quad_rule(n_nodes), data)
      fit_ses <- se_from_vcov(fit_vcov)
    } else {
      fit_vcov <- NA
      fit_ses <- NA
    }
  } else {
    stop("Invalid family specified. Please see documentations for available families.")
  }
  
  # prepare object for returning
  out <- list(
    family = family,
    fit = fit,
    fit_ses = fit_ses
  )
  
  return(out)
}