#' Estimate global shape and rate hyperparameters in BART-like fashion
#'
#' Uses a similar logic to that in Lamprinakou et. al (2023) to do a data
#' driven prior. It first tessellates the space (instead of partitioning via grid),
#' then treats the MLE of the intensity of each polygon divided by the number of learners
#' (or to the (1/number of learners) power if log_linear is TRUE)
#' as a 'sample'. The shape is set to the mean of the samples squared divided by the sample variance,
#' and the rate is set to the mean of the samples divided by the sample variance.
#'
#' @param boundary
#' @param data_points
#' @param N_LEARNERS
#' @param MESH_REF
#'
#' @return A numeric vector containing the shape and rate hyperparameters
#' @export
#'
#' @examples
bartStyleDataPrior = function(boundary, data_points, N_LEARNERS, MESH_REF, log_linear = FALSE) {
  data_prior_fit = makeVoronoiMeshes(boundary = boundary,
                                     N_LEARNERS = 1,
                                     MESH_REF = MESH_REF,
                                     max_area_ratio = 10) %>%
    makeLearners(data = data_points)

  shape_hp = 0
  rate_hp = 0
  if(log_linear) {
    emp_dens_samps = (data_prior_fit$node_counts[[1]])/(data_prior_fit$node_measures[[1]])
    emp_dens_candidates = emp_dens_samps^(1/N_LEARNERS)
    inv_empirical_var = 1/sd(emp_dens_candidates)^2
    shape_hp = (mean(emp_dens_candidates)^2)*(inv_empirical_var)
    rate_hp = (mean(emp_dens_candidates))*(inv_empirical_var)
  } else {
    emp_dens_samps = (data_prior_fit$node_counts[[1]])/(data_prior_fit$node_measures[[1]])
    emp_dens_candidates = emp_dens_samps/N_LEARNERS
    inv_empirical_var = 1/sd(emp_dens_candidates)^2
    shape_hp = (mean(emp_dens_candidates)^2)*(inv_empirical_var)
    rate_hp = (mean(emp_dens_candidates))*(inv_empirical_var)
  }

  return(c(shape_hp, rate_hp, inv_empirical_var))
}
