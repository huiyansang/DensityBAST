library(igraph)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(sf)
library(sp)
library(sn)
library(mgcv)
library(doParallel)
library(foreach)
library(doRNG)
library(mvtnorm)
library(TruncatedNormal)
library(scoringRules)
library(spatstat)
library(gridExtra)
library(furrr)
library(BART)

#' Fit DensityBAST model on linear point patterns
#'
#' @description
#' Fits the DensityBAST model for point patterns on linear networks using an unconditioned density approach.
#'
#' @param data_matrix An object representing the linear network point pattern (e.g., an lpp object).
#' @param hyperpars A list of hyperparameters for the model. If NULL, default hyperparameters are used:
#'   \itemize{
#'     \item N_LEARNERS = 16
#'     \item MESH_REF = 300
#'     \item MCMC = 2000
#'     \item BURNIN = 1000
#'     \item THIN = 1
#'     \item shape_hp = 1/16
#'     \item rate_hp = 1/16
#'     \item lambda_k_hp = 20
#'     \item max_clusters = 40
#'     \item inv_empirical_var = 1/16
#'     \item seed = NULL
#'     \item data_prior = TRUE
#'     \item normalized = FALSE
#'     \item hot_init = FALSE
#'     \item random_scan = FALSE
#'     \item detailed_output = TRUE
#'     \item subsetting = "thinning"
#'     \item thinning_ratio = 1
#'     \item parallel_cores = 16
#'   }
#'
#' @return A list containing the fitted DensityBAST model.
#' @export
DensityBASTFit_lpp=function(data_matrix,hyperpars=NULL) {
  if(is.null(hyperpars)){
    hyperpars = list(
      N_LEARNERS = 16,
      MESH_REF = 300,
      MCMC = 2000,
      BURNIN = 1000,
      THIN = 1,
      shape_hp = 1/16,
      rate_hp = 1/16,
      lambda_k_hp = 20,
      max_clusters = 40,
      inv_empirical_var = 1/16, # Just a default, this varies for each dataset!
      seed = NULL,
      data_prior = TRUE,
      normalized = FALSE,
      hot_init = FALSE,
      random_scan = FALSE,
      detailed_output = TRUE,
      subsetting = "thinning",
      thinning_ratio = 1,
      parallel_cores = 16
    )
  }

  test_linear_learners = makeLinearLearners(data_matrix = data_matrix,
                                            N_LEARNERS = hyperpars$N_LEARNERS,
                                            MESH_REF = hyperpars$MESH_REF)


  lin_dens_test_out = densityBASTunconditioned(
    test_linear_learners,
    MCMC = hyperpars$MCMC,
    BURNIN = hyperpars$BURNIN,
    THIN = hyperpars$THIN,
    shape_hp = hyperpars$shape_hp,
    rate_hp = hyperpars$rate_hp,
    lambda_k_hp = hyperpars$lambda_k_hp,
    max_clusters = hyperpars$max_clusters,
    inv_empirical_var = hyperpars$inv_empirical_var,
    seed = hyperpars$seed,
    data_prior = hyperpars$data_prior,
    normalized = hyperpars$normalized,
    hot_init = hyperpars$hot_init,
    random_scan = hyperpars$random_scan,
    subsetting = hyperpars$subsetting,
    thinning_ratio = hyperpars$thinning_ratio,
    detailed_output = hyperpars$detailed_output,
    parallel_cores = hyperpars$parallel_cores)

  return(bastion_fit= lin_dens_test_out)
}
#' Fit DensityBAST model on planar spatial point patterns
#'
#' @description
#' Fits the DensityBAST model for spatial point patterns within a given boundary. Returns functions
#' to estimate the density at arbitrary locations, including posterior mean and median estimates.
#'
#' @param boundary An sf polygon representing the spatial domain.
#' @param data_matrix A numeric matrix of point coordinates (columns = x, y).
#' @param hyperpars A list of hyperparameters. Defaults will be used if NULL.
#'
#' @return A list with the following functions:
#' \itemize{
#'   \item bastion_estFunc(coords) - returns density samples at the coordinates
#'   \item bastion_estFunc_postmean(coords) - returns posterior mean density at the coordinates
#'   \item bastion_estFunc_postmedian(coords) - returns posterior median density at the coordinates
#' }
#' @export
DensityBASTFit = function(boundary, data_matrix, hyperpars = NULL) {

  data_points = data_matrix %>%
    st_multipoint() %>%
    st_sfc %>%
    st_cast(to = "POINT")


  if(is.null(hyperpars)) {
    hyperpars = list(
      N_LEARNERS = 16,
      MESH_REF = 100,
      MCMC = 110,
      BURNIN = 10,
      THIN = 100,
      shape_hp = 1,
      rate_hp = 1,
      lambda_k_hp = 10,
      max_clusters = 20,
      inv_empirical_var = 1/16,
      subsetting = "thinning",
      thinning_ratio = 1/16,
      seed = NULL,
      data_prior = FALSE,
      normalized = FALSE,
      hot_init = FALSE,
      random_scan = FALSE,
      detailed_output = TRUE,
      parallel_cores = NULL
    )
  }

  if(hyperpars$data_prior) {
    prior_info = bartStyleDataPrior(boundary,
                                    data_points,
                                    N_LEARNERS = 1,
                                    MESH_REF = 1000)
    hyperpars$inv_empirical_var = prior_info[3]
  }

  densityOutput = with(hyperpars, {
    makeVoronoiMeshes(boundary = boundary,
                      N_LEARNERS = N_LEARNERS,
                      MESH_REF = MESH_REF,
                      max_area_ratio = 20) %>%
      makeLearners(data = data_points) %>%
      densityBASTunconditioned(learners = .,
                               MCMC = MCMC,
                               BURNIN = BURNIN,
                               THIN = THIN,
                               shape_hp = shape_hp,
                               rate_hp = rate_hp,
                               lambda_k_hp = lambda_k_hp,
                               max_clusters = max_clusters,
                               inv_empirical_var = inv_empirical_var,
                               seed = seed,
                               data_prior = data_prior,
                               normalized = normalized,
                               hot_init = hot_init,
                               random_scan = random_scan,
                               subsetting = subsetting,
                               thinning_ratio = thinning_ratio,
                               detailed_output = detailed_output,
                               parallel_cores = parallel_cores)
  })

  bastion_estFunc = function(coords) {
    coords %>%
      st_multipoint() %>%
      st_sfc %>%
      st_cast(to = "POINT") %>%
      st_sf %>%
      getDensitySamples(densityOutput, .)
  }

  bastion_estFunc_postmean = function(coords) {
    bastion_estFunc(coords) %>%
      colMeans()
  }

  bastion_estFunc_postmedian = function(coords) {
    bastion_estFunc(coords) %>%
      apply(2, median)
  }
  return(list(bastion_estFunc = bastion_estFunc,
              bastion_estFunc_postmean = bastion_estFunc_postmean,
              bastion_estFunc_postmedian = bastion_estFunc_postmedian))
}
#' Fit BART model for spatial point patterns
#'
#' @description
#' Fits a Bayesian Additive Regression Trees (BART) model for spatial point patterns using
#' a pseudo-likelihood approach with dummy points.
#'
#' @param boundary An sf polygon representing the spatial domain.
#' @param data_matrix A numeric matrix of point coordinates (columns = x, y).
#' @param hyperpars A list of hyperparameters (currently not fully utilized in this function).
#'
#' @return A list with the following function:
#' \itemize{
#'   \item bart_estFunc(coords) - returns estimated intensity at the coordinates
#' }
#' @export
DensityBartFit = function(boundary, data_matrix, hyperpars = NULL) {

  real_points = data.frame(x = data_matrix[,1],
                           y = data_matrix[,2])
  num_obs_points = nrow(data_matrix)
  # "We suggest taking delta (dummy intensity) = 4*n(X)/|D|"
  exp_num_dummy_points = 4*num_obs_points
  dummy_points = spatstat.random::rpoispp((exp_num_dummy_points)/sf::st_area(boundary),
                                          win=boundary)
  num_dummy_points = dummy_points$n
  delta_func = num_dummy_points/sf::st_area(boundary)
  total_points = num_dummy_points + num_obs_points

  sim_response = c(rep(0, num_dummy_points), rep(1, num_obs_points))
  sim_locations = data.frame(x = c(dummy_points$x, real_points$x),
                             y = c(dummy_points$y, real_points$y))

  bartOutput = pbart(x.train = sim_locations, y.train = sim_response, ndpost = 100)

  bart_estFunc = function(coords) {

    # Proxy for P(Y(u) = 1) = p_hat (posterior mean)
    # P(Y(u) = 1) = lambda(u)/(lambda(u) + delta(u))
    # -->
    # lambda(u) = -(delta(u)*p_hat)/(p_hat - 1)

    testOutput = predict(bartOutput, newdata = coords)

    test_pred = testOutput$prob.test.mean
    func_est = (delta_func*test_pred)/(1-test_pred)

    return(func_est)
  }

  return(list(bart_estFunc = bart_estFunc))

}


