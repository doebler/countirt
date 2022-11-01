summary.cirtfit <- function(x) {
  
  # TODO be able to account for turboEM output if accelerators are used
  
  # TODO change when i allow combination of constraints, cant be if else
  if (!is.null(x$model$fixed_alphas)) {
    con <- paste0("Slopes fixed to: ", 
                  paste(x$model$fixed_alphas, collapse = ", "))
  } else if (!is.null(x$model$fixed_log_disps)) {
    con <- paste0("Log dispersions fixed to: ", 
                  paste(x$model$fixed_log_disps, collapse = ", "))
  } else if (equal_alphas) {
    con <- "Slopes are constrained to be equal."
  } else if (equal_log_disps) {
    con <- "Log dispersions are constrained to be equal."
  }
  
  
  if (!is.null(x$model$i_cov_on)) {
    cov_type <- paste(x$model$i_cov_on, collapse = ", ")
  } else if (!is.null(x$model$p_covariates)) {
    cov_type <- "theta"
  } else {
    cov_type <- "none"
  }
  
  if (is.null(x$fit_ses)) {
    
  }
  # TODO make a data frame of the results so i can print that as the summary
  
  cat(paste0("Family: ", x$family, "\n"))
  cat(paste0("Constraints: ", con, "\n"))
  cat(paste0("Covariates: ", cov_type, "\n"))
  cat(paste0("Items: ", ncol(x$model$item_data), "\n"))
  cat(paste0("Observations: ", nrow(x$model$item_data), "\n"))
  cat("\n")
  cat("----------------------------------------------------------")
  cat("Model Parameters: \n")
  cat("----------------------------------------------------------")
  cat("Slopes: \n")
  
  cat("----------------------------------------------------------")
  cat("Intercepts: \n")
  
  if (x$family == "cmp") {
    cat("----------------------------------------------------------")
    cat("Log Dispersions: \n")
  }
  
}







