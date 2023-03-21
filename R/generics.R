#' Summary for cirtfit objects.
#' 
#' @param x A cirtfit object.
#' 
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @useDynLib countirt, .registration=TRUE
#' @export
summary.cirtfit <- function(x) {
  
  # TODO change whenever i allow new combination of constraints
  if (!is.null(x$model$fixed_alphas)) {
    con <- paste0("Slopes fixed to: ", 
                  paste(x$model$fixed_alphas, collapse = ", "))
  } else if (!is.null(x$model$fixed_log_disps)) {
    if (x$model$equal_alphas) {
      con <- paste0("Log dispersions fixed to: ", 
                    paste(x$model$fixed_log_disps, collapse = ", "), "\n",
                    "     Slopes are constrained to be equal.")
    } else {
      con <- paste0("Log dispersions fixed to: ", 
                    paste(x$model$fixed_log_disps, collapse = ", "))
    }
  } else if (x$model$equal_alphas) {
    con <- "Slopes are constrained to be equal."
  } else if (x$model$equal_log_disps) {
    con <- "Log dispersions are constrained to be equal."
  } else {
    # No constraints
    con <- "No constraints specified."
  }
  
  
  if (!is.null(x$model$i_cov_on)) {
    cov_type <- paste(x$model$i_cov_on, collapse = ", ")
  } else if (!is.null(x$model$p_covariates)) {
    cov_type <- "theta"
  } else {
    cov_type <- "none"
  }
  
  params <- data.frame(
    param = names(x$fit$params)
    )
  n_alphas <- sum(grepl("alpha", names(x$fit$params)))
  n_deltas <- sum(grepl("delta", names(x$fit$params)))
  n_logdisps <- sum(grepl("log_disp", names(x$fit$params)))
  params$Item <- c(rep("Global alpha", n_alphas), 
                   colnames(x$model$item_data), 
                   rep("Global log disp", n_logdisps))
  # TODO change if I allow more constraints
  if (n_alphas == ncol(x$model$item_data)) {
    # one alpha for each item
    params$Item[1:n_alphas] <- colnames(x$model$item_data)
  }
  if (n_logdisps == ncol(x$model$item_data)) {
    # one alpha for each item
    params$Item[(n_alphas+n_deltas+1):(n_alphas+n_deltas+n_logdisps)] <- colnames(x$model$item_data)
  }
  params$Estimate <- x$fit$params
  
  if (is.null(x$fit_ses)){
    params$`Est.Error` <- NA
    params$`z value` <- NA
    params$`p value` <- NA
    params$Significance <- NA
    warning("No standard errors, CIs, or z-tests were computed yet. Please use add_inference to add them if desired.")
  } else {
    params$`Est.Error` <- x$fit_ses$se
    params$`z value` <- x$fit_ses$z_value
    params$`p value` <- x$fit_ses$p_value
    params$Significance <- ifelse(
      x$fit_ses$p_value < 0.001, "***",
      ifelse(
        x$fit_ses$p_value < 0.01, "**",
        ifelse(
          x$fit_ses$p_value < 0.05, "*",
          ifelse(
            x$fit_ses$p_value < 0.1, ".",
            ""
          )
        )
      )
    )
  }
  
  params[,3:6] <- apply(params[,3:6], 2, round, digits = 3)
  
  if (!x$fit$conv) {
    warning("Model failed to converge.")
    conv <- "Model failed to converge."
  } else{
    conv <- paste0("Model converged after ", x$fit$iter, " EM iterations.")
  }
  
    cat(paste0(" Family: ", x$family, "\n"),
        paste0("Constraints: ", con, "\n"),
        paste0("Covariates on: ", cov_type, "\n"),
        paste0("Marginal log likelihood: ", round(x$fit$marg_ll[length(x$fit$marg_ll)],2), "\n"),
        paste0("Convergence: ", conv, "\n"),
        paste0("Items: ", ncol(x$model$item_data), "\n"),
        paste0("Observations: ", nrow(x$model$item_data), "\n"),
        "----------------------------------------------------------",
        "\n",
        "Model Parameters: \n",
        "----------------------------------------------------------",
        "\n")
    if (is.null(x$model$fixed_alphas)) {
      cat("Slopes: \n",
          "----------------------------------------------------------",
          "\n")
      print(params[grepl("alpha", params$param),-1])
      cat("\n",
          "----------------------------------------------------------",
          "\n")
    }
    cat("Intercepts: \n",
        "----------------------------------------------------------",
        "\n")
    print(params[grepl("delta", params$param),-1])
    if ((x$family == "cmp") & is.null((x$model$fixed_log_disps))) {
      cat("\n",
          "----------------------------------------------------------",
          "\n",
          "Log Dispersions: \n",
          "----------------------------------------------------------",
          "\n")
      print(params[grepl("log_disp", params$param),-1])
    }
    cat("\n")
    cat("----------------------------------------------------------")
    cat("\n")
    cat("Significance: '***' < .001, '**' < .01, '*' < .05, '.' < .1")
    cat("\n")
    cat("----------------------------------------------------------")
  invisible(x)
}

#' Compare two nested cirtfit object.
#' 
#' @param object A cirtfit object.
#' 
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling
#' @useDynLib countirt, .registration=TRUE
#' @export
anova.cirtfit <- function(object, ...) {
  
}





