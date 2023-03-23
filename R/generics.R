#' Summary for cirtfit objects.
#' 
#' @param x A cirtfit object.
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
  invisible(x)
}

#' Compare two nested cirtfit objects
#' 
#' Compare two nested cirtfit objects in using a likelihood ratio test, the AIC and the BIC.
#' 
#' @param object A cirtfit object.
#' @param ... Additional cirtfit objects.
#' 
#' @details The degrees of freedom used for the AIC and BIC are the number of model parameters. For the LR test, the degrees of freedom are the difference in the number of model parameters between the two models.
#' 
#' @useDynLib countirt, .registration=TRUE
#' @export
anova.cirtfit <- function(object, ...) {
  
  # TODO this only works for comparing two models at the moment
  # TODO make this work with explanatory models
  
  # make the other inputted models accessible
  models <- list(...)
  
  # check which one is the smaller model
  if (length(object$fit$params) < length(models[[1]]$fit$params)) {
    fit_small <- object
    fit_big <- models[[1]]
  } else {
    fit_small <- models[[1]]
    fit_big <- object
  }
  
  # LR test
  lr <- (-2)*(fit_small$fit$marg_ll[length(fit_small$fit$marg_ll)] - 
                fit_big$fit$marg_ll[length(fit_big$fit$marg_ll)])
  df <- length(fit_big$fit$params) - length(fit_small$fit$params)
  p_value <- 1 - pchisq(lr, df)
  
  # AIC 
  aic_big <- 2*length(fit_big$fit$params) - 2*fit_big$fit$marg_ll[length(fit_big$fit$marg_ll)]
  aic_small <- 2*length(fit_small$fit$params) - 2*fit_small$fit$marg_ll[length(fit_small$fit$marg_ll)]
  
  # BIC 
  n <- nrow(fit_small$model$item_data)
  bic_big <- length(fit_big$fit$params)*log(n) - 2*fit_big$fit$marg_ll[length(fit_big$fit$marg_ll)]
  bic_small <- length(fit_small$fit$params)*log(n) - 2*fit_small$fit$marg_ll[length(fit_small$fit$marg_ll)]
  
  
  # extract constraints for 
  # TODO adapt if i ever allow more constraints
  # TODO adapt to explanatory models
  if (!is.null(fit_big$model$fixed_alphas)) {
    con_big <- paste0("Slopes fixed to: ", 
                  paste(fit_big$model$fixed_alphas, collapse = ", "))
  } else if (!is.null(fit_big$model$fixed_log_disps)) {
    if (fit_big$model$equal_alphas) {
      con_big <- paste0("Log dispersions fixed to: ", 
                    paste(fit_big$model$fixed_log_disps, collapse = ", "), "\n",
                    "     Slopes are constrained to be equal.")
    } else {
      con_big <- paste0("Log dispersions fixed to: ", 
                    paste(fit_big$model$fixed_log_disps, collapse = ", "))
    }
  } else if (fit_big$model$equal_alphas) {
    con_big <- "Slopes are constrained to be equal."
  } else if (fit_big$model$equal_log_disps) {
    con_big <- "Log dispersions are constrained to be equal."
  } else {
    # No constraints
    con_big <- "No constraints specified."
  }
  
  if (!is.null(fit_small$model$fixed_alphas)) {
    con_small <- paste0("Slopes fixed to: ", 
                      paste(fit_small$model$fixed_alphas, collapse = ", "))
  } else if (!is.null(fit_small$model$fixed_log_disps)) {
    if (fit_small$model$equal_alphas) {
      con_small <- paste0("Log dispersions fixed to: ", 
                        paste(fit_small$model$fixed_log_disps, collapse = ", "), "\n",
                        "     Slopes are constrained to be equal.")
    } else {
      con_small <- paste0("Log dispersions fixed to: ", 
                        paste(fit_small$model$fixed_log_disps, collapse = ", "))
    }
  } else if (fit_small$model$equal_alphas) {
    con_small <- "Slopes are constrained to be equal."
  } else if (fit_small$model$equal_log_disps) {
    con_small <- "Log dispersions are constrained to be equal."
  } else {
    # No constraints
    con_small <- "No constraints specified."
  }
  
  # TODO im Output die modellbeschreibung einfuegen
  out <- data.frame(
    Model = c("Constr. Model", "Model"),
    Params = c(length(fit_small$fit$params), length(fit_big$fit$params)),
    MLL = c(round(fit_small$fit$marg_ll[length(fit_small$fit$marg_ll)],2), 
            round(fit_big$fit$marg_ll[length(fit_big$fit$marg_ll)],2)),
    AIC = c(aic_small, aic_big),
    BIC = c(bic_small, aic_big),
    Df = c("", df),
    ChiSqrd = c("", round(lr,2)),
    p = c("", round(p_value, 5)),
    Signif = c("",ifelse(
      p_value < 0.001, "***",
      ifelse(
        p_value < 0.01, "**",
        ifelse(
          p_value < 0.05, "*",
          ifelse(
            p_value < 0.1, ".",
            ""
          )
        )
      )
    )
  ))
  cat("Model Comparison")
  cat("\n")
  cat("-------------------------------------------------------------------------")
  cat("\n")
  cat(paste0("Family: ", object$family, "\n"))
  cat(paste0("Items: ", ncol(object$model$item_data), "\n"))
  cat("-------------------------------------------------------------------------")
  cat("\n")
  cat(paste0(paste0("Model: ", con_big, "\n")))
  cat(paste0(paste0("Constr. Model: ", con_small, "\n")))
  cat("-------------------------------------------------------------------------")
  cat("\n")
  print(out)
  cat("\n")
  cat("-------------------------------------------------------------------------")
  cat("\n")
  cat("Significance: '***' < .001, '**' < .01, '*' < .05, '.' < .1")
  
  invisible(object)
}





