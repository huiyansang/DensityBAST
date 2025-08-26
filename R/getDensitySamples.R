# Get (number of grid points) x (number of posterior samples) matrix for CRPS
# Each row corresponds to a grid point and each column is a posterior sample for that grid point

#' Get matrix of posterior samples for a given grid
#'
#' @param densityBASTOutput The return object from the densityBAST function
#' @param gridPOINTs An sf object containing each grid point in its own row in the geometry column
#' @param intensity Optional; defaults to TRUE. A boolean value specifying whether to return the intensity function values (TRUE) or normalize the values so that each posterior sample is a density (FALSE)
#'
#' @return A matrix with dimension (number of grid points) x (number of posterior samples); each row corresponds to a grid point and each column is a posterior sample for that grid point.
#' @export
#'
#' @examples
getDensitySamples = function(densityBASTOutput, gridPOINTs, intensity = TRUE, log_linear = FALSE) {
  lambda_output = densityBASTOutput$lambda_output
  membership_output = densityBASTOutput$membership_output
  meshes = densityBASTOutput$meshes
  n_post_samples = length(membership_output)
  n_grid_points = nrow(gridPOINTs)
  n_polygons = length(membership_output[[1]][[1]])
  n_learners = length(meshes)

  post_grid_clusts = matrix(nrow = n_grid_points, ncol = n_learners)
  density_samples = matrix(rep(0, n_grid_points*n_post_samples), nrow = n_post_samples, ncol = n_grid_points)
  grid_polygon_membership = vector("list", n_learners)

  for(learner in 1:n_learners) {
    learner_polygons = meshes[[learner]] %>% sf::st_sf(.)
    grid_polygon_membership[[learner]] = sf::st_nearest_feature(gridPOINTs, learner_polygons)
  }

  grid_poly_df = as.data.frame(grid_polygon_membership, col.names = 1:n_learners)

  if(log_linear) {
    log_density_samples = density_samples
    for(post_sample in 1:n_post_samples) {
      for(learner in 1:n_learners) {
        poly_constants = get_polygon_constant(membership_output[[post_sample]][[learner]], lambda_output[[post_sample]][[learner]])
        log_density_samples[post_sample, ] = log_density_samples[post_sample, ] + log(poly_constants[grid_poly_df[, learner]])
      }
      density_samples[post_sample, ] = exp(log_density_samples[post_sample, ])
    }
  } else {
    for(post_sample in 1:n_post_samples) {
      for(learner in 1:n_learners) {
        poly_constants = get_polygon_constant(membership_output[[post_sample]][[learner]], lambda_output[[post_sample]][[learner]])
        density_samples[post_sample, ] = density_samples[post_sample, ] + poly_constants[grid_poly_df[, learner]]
      }
    }
  }


  if(intensity == TRUE) {
    return(density_samples)
  } else {
    return(t(apply(density_samples, 1, function(x) x/sum(x))))
  }

}
