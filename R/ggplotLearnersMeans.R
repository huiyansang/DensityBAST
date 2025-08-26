#' Generate a list of ggplot's for each weak learner
#'
#' @param densityOutput The output of the densityBAST function
#' @param data Optional; if a single
#'
#' @return A list of length N_LEARNERS where each element is a ggplot of the posterior mean (considering only the samples in the corresponding learner)
#' @export
#'
#' @examples
ggplotLearnersMeans = function(densityOutput,
                               full_data,
                               partitions = NULL,
                               plotting_limit = 12) {

  N_LEARNERS = length(densityOutput$meshes)
  if(is.null(plotting_limit)) {
    plotting_limit = N_LEARNERS
  }

  ggplot_list = vector("list", plotting_limit)
  plotting_sf = vector("list", plotting_limit)

  for(learner in 1:plotting_limit) {
    posterior_means = rep(0, length(densityOutput$membership_output[[1]][[learner]]))

    for(sample in 1:length(densityOutput$lambda_output)) {
      # create a lookup table
      value_map = data.frame(cluster = seq(1, length(densityOutput$lambda_output[[sample]][[learner]])),
                             lambda = densityOutput$lambda_output[[sample]][[learner]])
      posterior_means = posterior_means + value_map$lambda[match(densityOutput$membership_output[[sample]][[learner]], value_map$cluster)]
    }
    posterior_means = posterior_means/length(densityOutput$lambda_output)

    plotting_sf[[learner]] = densityOutput$meshes[[learner]] %>% st_sf
    plotting_sf[[learner]]$post_mean = posterior_means
    if(is.null(partitions)) {
      ggplot_list[[learner]] = ggplot() +
        geom_sf(data = plotting_sf[[learner]], aes(fill = post_mean)) +
        scale_fill_viridis_c(limits = c(0, NA)) +
        geom_sf(data = full_data, size = 0.1)
    } else {
      ggplot_list[[learner]] = ggplot() +
        geom_sf(data = plotting_sf[[learner]], aes(fill = post_mean)) +
        scale_fill_viridis_c(limits = c(0, NA)) +
        geom_sf(data = full_data[which(partitions[[learner]])], size = 0.1)
    }
  }
  return(ggplot_list)
}
