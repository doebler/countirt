#' Model fitting function for count data IRT models
#' 
#' @param data data matrix, each column must correspond to one item, each row to an observational unit (usually a person)
#' @param family a string indicating the count data family, can be either "cmp" or "poisson"
#' @param n_nodes integer, number of quadrature nodes, defaults to 121, for the 2CMP model, no less than 100 quadrature nodes are recommended
#' @param stand_errors boolean, indicates whether standard errors for model parameters should be estimated, defaults to FALSE
#' @param constraints a list, indicating the constraints for the model, note that at the current moment, the possible constraints are mutually exclusive (this functionality will be extended in the future). 
#' @param control a list, providing control parameters for the estimation
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom fastGHQuad gaussHermiteData
#' @importFrom nleqslv nleqslv
#' @importFrom rootSolve gradient
#' @useDynLib countirt, .registration=TRUE
#' @export
cirt <- function(data, 
                 family,
                 n_nodes = 121,
                 stand_errors = FALSE,
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
  
  if (family == "cmp") {
    start_values <- get_start_values(
      data = data, init_disp_one = control$init_disp_one
    )
    
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
    
    }
  } else if (family == "poisson") {
    
  }
  
  
}