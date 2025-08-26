#' Perform BASTION Density Estimation with Unconditional Weak Learners
#'
#' @inheritParams densityBAST
#' @param multinomial Optional; if TRUE then 'data' will be properly partitioned rather than bootstrapped for each weak learner
#' @param parallel_cores Optional; the number of cores to register for a parallel cluster or NULL to not run in parallel
#'
#' @importFrom doRNG %dorng%
#'
#' @return
#' @export
#'
#' @examples
densityBASTunconditioned = function(learners,
                                    MCMC,
                                    BURNIN,
                                    THIN,
                                    shape_hp = NULL,
                                    rate_hp,
                                    lambda_k_hp,
                                    max_clusters,
                                    inv_empirical_var = 1,
                                    seed = NULL,
                                    data_prior = TRUE,
                                    normalized = TRUE,
                                    hot_init = FALSE,
                                    random_scan = FALSE,
                                    subsetting = c("multinomial", "dirichlet", "proper_partition", "thinning"),
                                    thinning_ratio = NULL,
                                    parallel_cores = NULL,
                                    detailed_output = FALSE,
                                    graph_detailed_output = FALSE) {


  if(!is.null(seed)) {
    set.seed(seed)
  }

  N_LEARNERS = length(learners$full_graphs)

  # Get information from the learners
  data = learners$data

  N = length(data)
  n_post_samples = (MCMC - BURNIN)#/THIN
  unconditionedOutputs = vector("list", N_LEARNERS)
  learner_data_subsets = vector("list", N_LEARNERS)
  partition = vector("list", N_LEARNERS)
  custom_mesh = vector("list", N_LEARNERS)
  lambda_correction = NULL

  if(N_LEARNERS != 1) {
    if(subsetting == "multinomial") {
      for(learner in 1:N_LEARNERS) {
        # Frequentist 'multinomial' bootstrap
        # In each learner, points sampled with equal prob
        included_indices = sample(1:N, size = N, replace = TRUE)
        partition[[learner]] = rep(0, N)
        partition[[learner]][included_indices] = learner
        learner_data_subsets[[learner]] = subset(data, (included_indices %in% 1:N))
      }
      lambda_correction = N_LEARNERS
    }

    if(subsetting == "dirichlet") {
      # Bayesian 'bootstrap' partition of the data
      # In each learner, points sampled with different prob
      for(learner in 1:N_LEARNERS) {
        included_indices = sample(1:N,
                                  size = N,
                                  replace = TRUE,
                                  prob = igraph::sample_dirichlet(1, rep(1, N)))
        partition[[learner]] = rep(0, N)
        partition[[learner]][included_indices] = learner
        learner_data_subsets[[learner]] = subset(data, (included_indices %in% 1:N))
      }
      lambda_correction = N_LEARNERS
    }

    if(subsetting == "proper_partition") {
      proper_partition = sample(1:N_LEARNERS,
                                size = N,
                                replace = TRUE)
      for(learner in 1:N_LEARNERS) {
        included_indices = which(proper_partition == learner)
        partition[[learner]] = proper_partition
        learner_data_subsets[[learner]] = subset(data, (included_indices %in% 1:N))
      }
      # Every point is only included once so no correction needed for additive form
      lambda_correction = 1
    }

    if(subsetting == "thinning") {
      if(is.null(thinning_ratio)) {
        thinning_ratio = 1/N_LEARNERS
      }
      else {
        if(thinning_ratio > 1) {
          warning("Thinning ratio should be <= 1; typically it is 1/N_LEARNERS.\n")
        }
      }
      for(learner in 1:N_LEARNERS) {
        included_indices = sample(c(FALSE, TRUE),
                                  size = N,
                                  replace = TRUE,
                                  prob = c(1-thinning_ratio, thinning_ratio)) %>% which()
        partition[[learner]] = rep(0, N)
        partition[[learner]][included_indices] = learner
        learner_data_subsets[[learner]] = subset(data, (included_indices %in% 1:N))
      }
      lambda_correction = N_LEARNERS*thinning_ratio
    }

  } else {
    included_indices = 1:N
    partition[[1]] = rep(0, N)
    partition[[1]][included_indices] = 1
    learner_data_subsets[[1]] = subset(data, (included_indices %in% 1:N))
    lambda_correction = 1
  }

  if(is.null(lambda_correction)) {
    stop("For more than 1 weak learner a valid subsetting method must be defined.\n
         Valid methods are 'multinomial', 'dirichlet', and 'proper_partition'.\n")
  }

  learners_list = vector("list", length = N_LEARNERS)

  if(thinning_ratio == 1 & subsetting == "thinning") {
    for(learner in 1:N_LEARNERS) {
      #node_contains = sf::st_contains(learners$meshes[[learner]], learner_data_subsets[[learner]])
      learners_list[[learner]] = list(
        full_graphs = learners$full_graphs[learner],
        node_measures = learners$node_measures[learner],
        node_contains = learners$node_contains[learner],
        node_counts = learners$node_counts[learner],
        data_membership = learners$data_membership[learner],
        graph_sizes = learners$graph_sizes[learner],
        data = learner_data_subsets[[learner]],
        meshes = learners$meshes[learner]
      )
    }
  } else {
    for(learner in 1:N_LEARNERS) {
      node_contains = sf::st_contains(learners$meshes[[learner]], learner_data_subsets[[learner]])
      learners_list[[learner]] = list(
        full_graphs = learners$full_graphs[learner],
        node_measures = learners$node_measures[learner],
        node_contains = list(node_contains),
        node_counts = list(sapply(node_contains, length)),
        data_membership = list(get_polygon_id(node_contains)),
        graph_sizes = learners$graph_sizes[learner],
        data = learner_data_subsets[[learner]],
        meshes = learners$meshes[learner]
      )
    }
  }



  if(is.null(parallel_cores)) {
    for(learner in 1:N_LEARNERS) {
      unconditionedOutputs[[learner]] = densityBAST(learners_list[[learner]],
                                                    MCMC,
                                                    BURNIN,
                                                    THIN,
                                                    shape_hp,
                                                    rate_hp,
                                                    lambda_k_hp,
                                                    max_clusters,
                                                    inv_empirical_var,
                                                    seed = NULL,
                                                    data_prior = data_prior,
                                                    normalized = normalized,
                                                    random_scan = random_scan,
                                                    hot_init = hot_init,
                                                    detailed_output = detailed_output,
                                                    graph_detailed_output = graph_detailed_output)
    }
  } else {
    BASTIONcluster = parallel::makeCluster(parallel_cores)
    doParallel::registerDoParallel(BASTIONcluster)

    unconditionedOutputs = foreach::foreach(learner = 1:N_LEARNERS, .packages = "BASTION") %dorng% {
      densityBAST(learners_list[[learner]],
                  MCMC,
                  BURNIN,
                  THIN,
                  shape_hp,
                  rate_hp,
                  lambda_k_hp,
                  max_clusters,
                  inv_empirical_var,
                  seed = NULL,
                  data_prior = data_prior,
                  normalized = normalized,
                  random_scan = random_scan,
                  hot_init = hot_init,
                  detailed_output = detailed_output,
                  graph_detailed_output = graph_detailed_output)
    }

    parallel::stopCluster(BASTIONcluster)
  }


  if(detailed_output & graph_detailed_output) {
    unconditionedCombinedOutput = list("lambda_output" = vector("list", n_post_samples),
                                       "membership_output" = vector("list", n_post_samples),
                                       "meshes" = vector("list", N_LEARNERS),
                                       "graph_output" = vector("list", n_post_samples),
                                       "acceptance_output" = vector("list", n_post_samples),
                                       "log_marginal_output" = vector("list", n_post_samples),
                                       "point_inclusion_output" = vector("list", n_post_samples),
                                       "partitions_output" = vector("list", N_LEARNERS))
    unconditionedCombinedOutput$partitions_output = partition
  } else if(detailed_output & !graph_detailed_output) {
    unconditionedCombinedOutput = list("lambda_output" = vector("list", n_post_samples),
                                       "membership_output" = vector("list", n_post_samples),
                                       "meshes" = vector("list", N_LEARNERS),
                                       "acceptance_output" = vector("list", n_post_samples),
                                       "log_marginal_output" = vector("list", n_post_samples),
                                       "point_inclusion_output" = vector("list", n_post_samples),
                                       "partitions_output" = vector("list", N_LEARNERS))
    unconditionedCombinedOutput$partitions_output = partition
  } else {
    unconditionedCombinedOutput = list("lambda_output" = vector("list", n_post_samples),
                                       "membership_output" = vector("list", n_post_samples),
                                       "meshes" = vector("list", N_LEARNERS))
  }


  for(learner in 1:N_LEARNERS) {
    unconditionedCombinedOutput$meshes[[learner]] = unconditionedOutputs[[learner]]$meshes[[1]]
  }

  for(post_sample in 1:n_post_samples) {
    for(learner in 1:N_LEARNERS) {
      # Divide by lambda_correction to correct intensity for each learner having (lambda_correction) points
      unconditionedCombinedOutput$lambda_output[[post_sample]][[learner]] = unconditionedOutputs[[learner]]$lambda_output[[post_sample]][[1]]/lambda_correction
      unconditionedCombinedOutput$membership_output[[post_sample]][[learner]] = unconditionedOutputs[[learner]]$membership_output[[post_sample]][[1]]
      if(detailed_output) {
        #
        unconditionedCombinedOutput$acceptance_output[[post_sample]][[learner]] = unconditionedOutputs[[learner]]$acceptance_output[[post_sample]][[1]]
        unconditionedCombinedOutput$log_marginal_output[[post_sample]][[learner]] = unconditionedOutputs[[learner]]$log_marginal_output[[post_sample]][[1]]
        unconditionedCombinedOutput$point_inclusion_output[[post_sample]][[learner]] = unconditionedOutputs[[learner]]$point_inclusion_output[[post_sample]][[1]]
      }
      if(graph_detailed_output) {
        unconditionedCombinedOutput$graph_output[[post_sample]][[learner]] = unconditionedOutputs[[learner]]$graph_output[[post_sample]][[1]]
      }
    }
  }

  return(unconditionedCombinedOutput)

}
