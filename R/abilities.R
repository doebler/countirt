#' Extract ability parameter estimates from a cirtfit object
#' 
#' Extract ability parameter estimates from a cirtfit object. There are different possibiliries for ability parameter estimation; currently, the only available option is EAP estimates.
#' 
#' @param x A cirtfit object.
#' @param type Defaults to "EAP" which is the only currently available option.
#' 
#' @useDynLib countirt, .registration=TRUE
#' @export
abilities <- function(x, type = "EAP") {
  
  if (type == "EAP") {
    if (x$family == "cmp") {
      out <- get_ability_params_pp(
        x$model$item_data, 
        x$fit$params, 
        x$control$n_nodes)
    } else if (x$family == "poisson") {
      out <- get_ability_params_poisson_pp(
        x$model$item_data, 
        x$fit$params, 
        x$control$n_nodes
      )
    }
  } else {
    stop("Invalid type specified.")
  }
  
  return(out)
}