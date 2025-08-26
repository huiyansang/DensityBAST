#' Evaluate performance of densityBAST against ground truth
#'
#' This function prints out Root Integrated Square Error, Average Absolute Error for
#' both the posterior mean and median as point estimates for the estimated intensity
#' function; as well as with just the posterior median for the estimated density function
#'
#' @param density_output The output of a densityBAST function
#' @param grid_df A 3 column data frame containing a grid of 'X', 'Y', and 'true_intens' true intensity values at those coordinates
#' @param log_linear Passed into getDensitySamples; whether to interpret the output as coming from a log-linear model
#' @param printout Whether or not to print results to the console
#'
#' @return A list of containing plot_list, error_list, and est_intense_list
#' @export
#'
#'
#' @examples
summarizeDensityFit = function(density_output,
                               grid_df,
                               log_linear = FALSE,
                               printout = TRUE) {

  if(any(colnames(grid_df) != c("X", "Y", "true_intens"))) {
    stop("grid_df colnames must be 'X', 'Y', 'true_intens'!")
  }

  true_intens = grid_df[,3]

  grid_points = (grid_df[,1:2]) %>%
    as.matrix() %>%
    sf::st_multipoint() %>%
    sf::st_sfc() %>%
    sf::st_cast(., to = "POINT") %>%
    sf::st_sf()
  BASTION_intensity_samples = getDensitySamples(density_output,
                                               grid_points,
                                               intensity = TRUE,
                                               log_linear = log_linear)
  BASTION_density_samples = t(apply(BASTION_intensity_samples, 1,
                                         function(x) x/sum(x)))

  BASTION_mean_intense = colMeans(BASTION_intensity_samples)
  BASTION_mean_sqerr = (BASTION_mean_intense - true_intens)^2
  BASTION_mean_RISE = sqrt(mean(BASTION_mean_sqerr))
  BASTION_mean_AAE = mean(abs(BASTION_mean_intense - true_intens))


  BASTION_median_intense = apply(BASTION_intensity_samples, 2, median)
  BASTION_median_sqerr = (BASTION_median_intense - true_intens)^2
  BASTION_median_RISE = sqrt(mean(BASTION_median_sqerr))
  BASTION_median_AAE = mean(abs(BASTION_median_intense - true_intens))


  true_density = true_intens/sum(true_intens)
  BASTION_density = apply(BASTION_density_samples, 2, median)
  BASTION_density_sqerr = (BASTION_density - true_density)^2
  BASTION_density_RISE = sqrt(mean(BASTION_density_sqerr))
  BASTION_density_AAE = mean(abs(BASTION_density - true_density))

  # Plot the BASTION intensity estimate (posterior mean)
  post_mean_plot = ggplot2::ggplot(cbind(grid_df, BASTION_mean_intense)) +
    ggplot2::geom_raster(ggplot2::aes(x = X, y = Y, fill = BASTION_mean_intense)) +
    ggplot2::scale_fill_viridis_c() #+
  #geom_point(data = data_matrix %>% as.data.frame(),
  #           aes(x = V1, y = V2), size = 0.05)
  # Plot the Squared Error
  post_mean_sqerr_plot = ggplot2::ggplot(cbind(grid_df, BASTION_mean_sqerr)) +
    ggplot2::geom_raster(ggplot2::aes(x = X, y = Y, fill = BASTION_mean_sqerr)) +
    ggplot2::scale_fill_continuous() #+
  #geom_point(data = data_matrix %>% as.data.frame(),
  #           aes(x = V1, y = V2), size = 0.05)

  # Plot the BASTION intensity estimate (posterior median)
  post_median_plot = ggplot2::ggplot(cbind(grid_df, BASTION_median_intense)) +
    ggplot2::geom_raster(ggplot2::aes(x = X, y = Y, fill = BASTION_median_intense)) +
    ggplot2::scale_fill_viridis_c() #+
  #geom_point(data = data_matrix %>% as.data.frame(),
  #           aes(x = V1, y = V2), size = 0.05)
  # Plot the Squared Error
  post_median_sqerr_plot = ggplot2::ggplot(cbind(grid_df, BASTION_median_sqerr)) +
    ggplot2::geom_raster(ggplot2::aes(x = X, y = Y, fill = BASTION_median_sqerr)) +
    ggplot2::scale_fill_continuous() #+
  #geom_point(data = data_matrix %>% as.data.frame(),
  #           aes(x = V1, y = V2), size = 0.05)


  if(printout) {
    cat("The Root Integrated Squared Error of the BASTION predicted
    intensity function using posterior mean is:", BASTION_mean_RISE, "\n")
    cat("The Average Absolute Error of the BASTION predicted
    intensity function using posterior mean is:", BASTION_mean_AAE, "\n")

    cat("The Root Integrated Squared Error of the BASTION predicted
    intensity function using posterior median is:",
        BASTION_median_RISE, "\n")
    cat("The Average Absolute Error of the BASTION predicted
    intensity function using posterior median is:",
        BASTION_median_AAE, "\n")

    cat("The Root Integrated Squared Error of the BASTION predicted
    density function using posterior median is:",
        BASTION_density_RISE, "\n")
    cat("The Average Absolute Error of the BASTION predicted
    density function using posterior median is:",
        BASTION_density_AAE, "\n")
  }

  plot_list = list(post_mean_plot,
                   post_mean_sqerr_plot,
                   post_median_plot,
                   post_median_sqerr_plot)
  error_list = list(BASTION_mean_RISE,
                    BASTION_mean_AAE,
                    BASTION_median_RISE,
                    BASTION_median_AAE)
  est_intense_list = list(BASTION_mean_intense,
                          BASTION_median_intense)

  return(list(plot_list, error_list, est_intense_list))
}
