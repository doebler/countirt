#' Function to plot item response curves for cirtfit object
#' 
#' @param x A cirtfit object.
#' @param shading Defaults to TRUE. Determines whether model implied SD of conditional response distributions should be drawn into the plot as shaded band (+/- 1 SD).
#' @param grid Defaults to TRUE. Deermines whether all item curves should be shown in one panel or each in their own panel.
#' @param theta_grid Defaults to 300. The number of theta values for theta grid for x axis of plot. The finer the grid, the smoother the item curve.
#' 
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @export
item_curves <- function(x, shading = TRUE, grid = FALSE, theta_grid = 300) {
  # TODO add more arguments to customize plots more
  
  n_items <- ncol(x$model$item_data)
  item_names <- colnames(x$model$item_data)
  
  # check for constraints to see if we can just extract parameters or we
  # need to do something to them
  # TODO change whenever i allow new combination of constraints
  deltas <- x$fit$params[grepl("delta", names(x$fit$params))]
  if (!is.null(x$model$fixed_alphas)) {
    alphas <- x$model$fixed_alphas
    if (x$family == "cmp") {
      disps <- exp(x$fit$params[grepl("log_disp", names(x$fit$params))])
    }
  } else if (!is.null(x$model$fixed_log_disps)) {
    if (x$model$equal_alphas) {
      alphas <- rep(x$fit$params[grepl("alpha", names(x$fit$params))], n_items)
      disps <- exp(x$model$fixed_log_disps)
    } else {
      alphas <- x$fit$params[grepl("alpha", names(x$fit$params))]
      disps <- exp(x$model$fixed_log_disps)
    }
  } else if (x$model$equal_alphas) {
    alphas <- rep(x$fit$params[grepl("alpha", names(x$fit$params))], n_items)
    if (x$family == "cmp") {
      disps <- exp(x$fit$params[grepl("log_disp", names(x$fit$params))])
    }
  } else if (x$model$equal_log_disps) {
    alphas <- x$fit$params[grepl("alpha", names(x$fit$params))]
    disps <- rep(exp(x$fit$params[grepl("log_disp", names(x$fit$params))]), n_items)
  } else {
    # No constraints
    alphas <- x$fit$params[grepl("alpha", names(x$fit$params))]
    if (x$family == "cmp") {
      disps <- exp(x$fit$params[grepl("log_disp", names(x$fit$params))])
    }
  }
  
  # thetas for x axis
  item_params <- cbind(alphas, deltas)
  thetas <- seq(-3, 3, length.out = theta_grid)
  df <- matrix(NA, nrow = theta_grid, ncol = n_items)
  df <- apply(item_params, 1, function(x){exp(x[1]*thetas + x[2])})
  df <- as.data.frame(df)
  df$thetas <- thetas
  colnames(df) <- c(item_names, "thetas")

  # prepare the shading
  df_long <- pivot_longer(df, -thetas, names_to = "Item", values_to = "mus")
  if (shading) {
    if (x$family == "cmp") {
      df_disps <- data.frame(item_names, disps)
      disp_pos <- match(df_long$Item, df_disps[,1])
      df_long$disps <- sapply(disp_pos, function(x){df_disps[x,2]})
      df_long$vars <- countirt:::get_var_cmp(mu = df_long$mus, nu = df_long$disps)
    } else if (x$family == "poisson") {
      df_long$vars <- df_long$mus
    } else {
      stop("Invalid family specified.")
    }
    df_long$plus_sd <- df_long$mus + sqrt(df_long$vars)
    df_long$minus_sd <- df_long$mus - sqrt(df_long$vars)
  }
  
  # make the plot
  if (shading) {
    pl <- ggplot(data = df_long,
           aes(x = thetas, y = mus, 
               ymin = minus_sd, ymax = plus_sd,
               col = Item, fill = Item, linetype = Item)) +
      geom_line() +
      geom_ribbon(alpha=0.1) +
      xlab("Theta") + 
      ylab("Mu (Expected Counts)") +
      scale_color_viridis_d() +
      scale_fill_viridis_d() +
      theme_bw()
  } else {
    pl <- ggplot(data = df_long,
           aes(x = thetas, y = mus, 
               col = Item, linetype = Item)) +
      geom_line() +
      xlab("Theta") + 
      ylab("Mu (Expected Counts)") +
      scale_color_viridis_d() +
      theme_bw()
  }
  
  if (grid) {
    pl <- pl +
      facet_wrap(~ Item)
  }
  
  pl
}











