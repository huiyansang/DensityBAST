#' Perform BASTION Density Estimation
#'
#' @param data An object of class sfc_POINT; a list of the coordinates of observed data locations
#' @param boundary An object of class sfg POLYGON; the boundary which defines where the data can lie
#' @param N_LEARNERS An integer, the number of weak learners to use when training the model
#' @param MESH_REF An integer, the number of reference nodes in each weak learner's mesh. More reference nodes leads to a finer mesh
#' @param MCMC An integer, the number of iterations for the MCMC
#' @param BURNIN An integer, the number of initial iterations for the MCMC before recording output
#' @param THIN An integer, specifies the thinning interval
#' @param shape_hp Shape hyperparameter for the gamma prior. By default it is 1/2
#' @param rate_hp Rate hyperparameter for the gamma prior. Ignored if data_prior = TRUE
#' @param lambda_k_hp Hyperparameter for poisson prior on the number of clusters
#' @param max_clusters An integer, the maximum number of clusters for each weak learner
#' @param seed Optional, value to feed to R's set.seed function
#'
#' @importFrom magrittr %>%
#' @importFrom doRNG %dorng%
#'
#'
#' @return A list containing two elements:
#' \item{lambda_output}{The piecewise constants corresponding to each cluster in each learner}
#' \item{membership_output}{The cluster memberships of ecah data point in each learner}
#' @export
#'
#' @examples
densityBASTunblocked = function(learners,
                       MCMC,
                       BURNIN,
                       THIN,
                       shape_hp = 0.5,
                       rate_hp,
                       lambda_k_hp,
                       max_clusters,
                       seed = NULL,
                       data_prior = FALSE,
                       normalized = FALSE,
                       hot_init = FALSE,
                       random_scan = FALSE,
                       detailed_output = FALSE,
                       graph_detailed_output = FALSE,
                       parallel_cores = NULL) {

  if(!is.null(seed)) {
    set.seed(seed)
  }
  if(!is.null(parallel_cores)) {
    BASTIONcluster = parallel::makeCluster(parallel_cores)
    doParallel::registerDoParallel(BASTIONcluster)
  }

  # Number of learners
  N_LEARNERS = length(learners$meshes)
  # Number of observations
  N = length(learners$data)
  # Number of posterior samples
  n_post_samps = (MCMC-BURNIN)
  # Get information from the learners
  full_graphs = learners$full_graphs
  polygon_areas = learners$node_measures
  point_in_polygon = learners$node_contains
  points_per_polygon = learners$node_counts
  data_polygon_id = learners$data_membership
  n_polygons = learners$graph_sizes
  meshes = learners$meshes
  data = learners$data

  ## Setup section - storage initialization
  current_graphs = vector("list", N_LEARNERS)
  thinned_ppp_proposed = vector("list", N_LEARNERS)
  thinned_ppp_current = vector("list", N_LEARNERS)
  learner_components = vector("list", N_LEARNERS)
  lambda_values = vector("list", N_LEARNERS)
  polygon_constants = vector("list", N_LEARNERS)
  acceptance_prob = rep(0, N_LEARNERS)
  log_Lik = rep(0, N_LEARNERS)
  cluster_areas = vector("list", N_LEARNERS)
  cluster_ppp = vector("list", N_LEARNERS)

  # Output storage initialization
  lambda_output = vector("list", n_post_samps)
  membership_output = vector("list", n_post_samps)
  if(detailed_output) {
    acceptance_output = vector("list", n_post_samps)
    log_marginal_output = vector("list", n_post_samps)
    #thinned_ppp_output = vector("list", n_post_samps)
    point_inclusion_output = vector("list", n_post_samps)
  }
  if(graph_detailed_output) {
    graph_output = vector("list", n_post_samps)
  }
  # NOTE: not yet outputting posterior probabilities
  if(N_LEARNERS == 1) {
    do_thinning = FALSE
  } else {
    do_thinning = TRUE
  }

  ## Setup section - weak learner initialization
  for(learner in 1:N_LEARNERS) {

    # Scale the total area to N
    original_area = sum(polygon_areas[[learner]])
    if(normalized) {
      polygon_areas[[learner]] = polygon_areas[[learner]]*(N/original_area)
    }
    total_area = sum(polygon_areas[[learner]])
    current_graphs[[learner]] = igraph::mst(full_graphs[[learner]]) # Generate a minimum spanning tree
    points_per_polygon[[learner]] = sapply(point_in_polygon[[learner]], length) # Get the number of points in each polygon
    # Instead of starting with 1 cluster, start with the maximum number of clusters
    # Encourage mixing by doing a hyper before every split
    learner_components[[learner]] = igraph::components(current_graphs[[learner]]) # Membership & number of clusters in current graph
    k = learner_components[[learner]]$no

    if(hot_init) {
      for(splits in 1:max_clusters) {
        learner_components[[learner]] = igraph::components(current_graphs[[learner]])
        k = learner_components[[learner]]$no
        cluster_split_index = sample.int(k, 1, prob = (learner_components[[learner]]$csize - 1))
        proposed_graph = (graphHyper(full_graphs[[learner]],
                                    learner_components[[learner]]$membership))$graph %>%
                          graphBirth(.,
                                    learner_components[[learner]]$membership,
                                    cluster_split_index)
        current_graphs[[learner]] = proposed_graph$graph
      }
      learner_components[[learner]] = igraph::components(current_graphs[[learner]])
      k = learner_components[[learner]]$no
    }

    cluster_areas[[learner]] = rep(0, k)
    cluster_ppp[[learner]] = rep(0, k)

    if(do_thinning) {
      point_included = purrr::rbernoulli(N, p = (1/N_LEARNERS))
      thinned_ppp_current[[learner]] = data_polygon_id[[learner]][point_included] %>%
        factor(.,
               levels = 1:length(points_per_polygon[[learner]])) %>%
        table() %>%
        as.numeric()
    } else {
      thinned_ppp_current[[learner]] = points_per_polygon[[learner]]
    }


    for(clust in 1:(k)) {
      cluster_areas[[learner]][clust] = sum(polygon_areas[[learner]][learner_components[[learner]]$membership == clust])
      cluster_ppp[[learner]][clust] = sum(thinned_ppp_current[[learner]][learner_components[[learner]]$membership == clust])
    }

    lambda_values[[learner]] = generateLambdaVals(data_prior,
                                                  k,
                                                  shape_hp,
                                                  rate_hp,
                                                  cluster_ppp[[learner]],
                                                  cluster_areas[[learner]])

    polygon_constants[[learner]] = get_polygon_constant(learner_components[[learner]]$membership,
                                                        lambda_values[[learner]])
  }


  # Initialize progress bar
  pb <- progress::progress_bar$new(
    format = "  MCMC running [:bar] :percent eta: :eta",
    total = MCMC, clear = FALSE, width= 60)

  for(iter in 1:MCMC) {
    pb$tick()
    output_samp = iter-BURNIN
    for(rep in 1:THIN) {
      if(random_scan) {
        sampling_order = sample(1:N_LEARNERS, N_LEARNERS, replace = FALSE)
      } else {
        sampling_order = 1:N_LEARNERS
      }
      for(learner in sampling_order) {
        if(do_thinning == FALSE) {
          point_included = rep(TRUE, N)
        } else {
          # Determine thinning probability:
          rho_s = numeric(N)
          for(i in 1:N_LEARNERS) {
            rho_s = rho_s + polygon_constants[[i]][data_polygon_id[[i]]]
          }
          delta_s = (polygon_constants[[learner]][data_polygon_id[[learner]]])/rho_s
          # Thin the points according to the thinning probabilities
          point_included = purrr::rbernoulli(N, p = delta_s)
          if(sum(point_included) == 0) {
            warning(paste0("Zero'd out weak learner: ", learner, ", consider reducing N_LEARNERS\n"))
          }
          if(detailed_output & (output_samp > 0)) {
            point_inclusion_output[[output_samp]][[learner]] = point_included
          }
        }
        thinned_ppp_proposed[[learner]] = data_polygon_id[[learner]][point_included] %>%
          factor(.,
                 levels = 1:length(polygon_constants[[learner]])) %>%
          table() %>%
          as.numeric()
        # If do_thinning == FALSE, then this is just equal to points_per_polygon[[learner]]
        thinned_ppp_current[[learner]] = thinned_ppp_proposed[[learner]]

        # current number of clusters
        k = learner_components[[learner]]$no
        # define move probabilities
        if(k == 1) {
          rb = 0.9; rd = 0; rc = 0; rhy = 0.1
        } else if(k == min(max_clusters, n_polygons[[learner]])) {
          rb = 0; rd = 0.6; rc = 0.3; rhy = 0.1
        } else {
          rb = 0.3; rd = 0.3; rc = 0.3; rhy = 0.1
        }
        move = sample(4, 1, prob = c(rb, rd, rc, rhy))

        # birth move
        if(move == 1) {
          # Propose a split
          # -1 to prevent splitting cluster with only 1 member
          cluster_split_index = sample.int(k, 1, prob = (learner_components[[learner]]$csize - 1))
          proposed_graph = graphBirth(current_graphs[[learner]],
                                      learner_components[[learner]]$membership,
                                      cluster_split_index)
          # Compute the acceptance probability
          # compute log-proposal ratio
          if(k == min(max_clusters, N)-1) {
            rd_new = 0.6
          } else {
            rd_new = 0.3
          }
          log_Prop = log(rd_new) - log(rb)
          # compute the log-prior ratio on k; the number of clusters
          log_Prior = log(lambda_k_hp) - log(k + 1)
          # compute log-marginal likelihood ratio
          log_Lik[learner] = density_logFullCorrectedMarginal(current_graphs[[learner]],
                                                              proposed_graph$graph,
                                                              shape_hp,
                                                              rate_hp,
                                                              thinned_ppp_proposed[[learner]],
                                                              thinned_ppp_current[[learner]],
                                                              polygon_areas[[learner]],
                                                              data_prior = data_prior,
                                                              inf.warn = FALSE)
          acc_prob = min(0, log_Prop + log_Lik[learner] + log_Prior)
          acc_prob = exp(acc_prob)
          if(detailed_output) {
            acceptance_prob[learner] = acc_prob
          }
          if(is.nan(acc_prob)) {
            warning("Acceptance probability is NaN\n")
          }
          if(runif(1) < acc_prob){
            # accept
            current_graphs[[learner]] = proposed_graph$graph
            #thinned_ppp_current[[learner]] = thinned_ppp_proposed[[learner]]
          }

        }
        # death move
        if(move == 2) {
          # Propose a merge
          proposed_graph = graphDeath(current_graphs[[learner]],
                                      learner_components[[learner]]$membership,
                                      full_graphs[[learner]])
          # Compute the acceptance probability
          # compute log-proposal ratio
          if(k == 2) {
            rb_new = 0.9
          } else {
            rb_new = 0.3
          }
          log_Prop = -(log(rd) - log(rb_new))
          # compute log-prior ratio on k; the number of clusters
          log_Prior = log(k) - log(lambda_k_hp)
          # compute log-posterior ratio
          log_Lik[learner] = density_logFullCorrectedMarginal(current_graphs[[learner]],
                                                              proposed_graph$graph,
                                                              shape_hp,
                                                              rate_hp,
                                                              thinned_ppp_proposed[[learner]],
                                                              thinned_ppp_current[[learner]],
                                                              polygon_areas[[learner]],
                                                              data_prior = data_prior,
                                                              inf.warn = FALSE)

          acc_prob = min(0, log_Prop + log_Lik[learner] + log_Prior)
          acc_prob = exp(acc_prob)
          if(detailed_output) {
            acceptance_prob[learner] = acc_prob
          }
          if(is.nan(acc_prob)) {
            warning("Acceptance probability is NaN\n")
          }
          if(runif(1) < acc_prob) {
            # accept
            current_graphs[[learner]] = proposed_graph$graph
            #thinned_ppp_current[[learner]] = thinned_ppp_proposed[[learner]]
          }
        }
        # change move
        if(move == 3) {
          # Propose a merge then a split
          proposed_graph = graphChange(current_graphs[[learner]],
                                       learner_components[[learner]]$membership,
                                       full_graphs[[learner]])
          # compute log-posterior ratio
          log_Lik_d = density_logFullCorrectedMarginal(current_graphs[[learner]],
                                                       proposed_graph$intermediate_graph,
                                                       shape_hp,
                                                       rate_hp,
                                                       thinned_ppp_proposed[[learner]],
                                                       thinned_ppp_current[[learner]],
                                                       polygon_areas[[learner]],
                                                       data_prior = data_prior,
                                                       inf.warn = FALSE)
          log_Lik_b = density_logFullCorrectedMarginal(proposed_graph$intermediate_graph,
                                                       proposed_graph$graph,
                                                       shape_hp,
                                                       rate_hp,
                                                       thinned_ppp_proposed[[learner]],
                                                       thinned_ppp_current[[learner]],
                                                       polygon_areas[[learner]],
                                                       data_prior = data_prior,
                                                       inf.warn = FALSE)
          log_Lik[learner] = log_Lik_b + log_Lik_d
          acc_prob = min(0, log_Lik[learner])
          acc_prob = exp(acc_prob)
          if(detailed_output) {
            acceptance_prob[learner] = acc_prob
          }
          if(is.nan(acc_prob)) {
            warning("Acceptance probability is NaN\n")
          }
          if(runif(1) < acc_prob){
            # accept
            current_graphs[[learner]] = proposed_graph$graph
            #thinned_ppp_current[[learner]] = thinned_ppp_proposed[[learner]]
          }

        }
        # hyper move
        if(move == 4) {
          proposed_graph = graphHyper(full_graphs[[learner]], learner_components[[learner]]$membership)
          current_graphs[[learner]] = proposed_graph$graph
          if(detailed_output) {
            acceptance_prob[learner] = 1
          }
          log_Lik[learner] = 0
        }

        learner_components[[learner]] = igraph::components(current_graphs[[learner]])
        k = learner_components[[learner]]$no
        cluster_areas[[learner]] = rep(0, k)
        cluster_ppp[[learner]] = rep(0, k)

        for(clust in 1:(k)) {
          cluster_areas[[learner]][clust] = sum(polygon_areas[[learner]][learner_components[[learner]]$membership == clust])
          cluster_ppp[[learner]][clust] = sum(thinned_ppp_current[[learner]][learner_components[[learner]]$membership == clust])
        }
        lambda_values[[learner]] = generateLambdaVals(data_prior,
                                                      k,
                                                      shape_hp,
                                                      rate_hp,
                                                      cluster_ppp[[learner]],
                                                      cluster_areas[[learner]])

        polygon_constants[[learner]] = get_polygon_constant(learner_components[[learner]]$membership, lambda_values[[learner]])
      }
    }


    if(output_samp > 0) {
      for(learner in 1:N_LEARNERS) {
        if(normalized) {
          lambda_output[[output_samp]][[learner]] = lambda_values[[learner]]*(N/original_area)
        } else {
          lambda_output[[output_samp]][[learner]] = lambda_values[[learner]]
        }
        if(data_prior) {
          lambda_output[[output_samp]][[learner]] = 2*lambda_output[[output_samp]][[learner]]
        }
        membership_output[[output_samp]][[learner]] = learner_components[[learner]]$membership
        if(detailed_output) {
          acceptance_output[[output_samp]][[learner]] = acceptance_prob[learner]
          log_marginal_output[[output_samp]][[learner]] = log_Lik[learner]
          #thinned_ppp_output[[output_samp]][[learner]] = thinned_ppp_current[[learner]]
        }
        if(graph_detailed_output) {
          graph_output[[output_samp]][[learner]] = current_graphs[[learner]]
        }
      }
    }
  }

  output = list(lambda_output = lambda_output,
                membership_output = membership_output,
                meshes = meshes)

  if(detailed_output) {
    output$acceptance_output = acceptance_output
    output$log_marginal_output = log_marginal_output
    #output$thinned_ppp_output = thinned_ppp_output
    output$point_inclusion_output = point_inclusion_output
  }

  if(graph_detailed_output) {
    output$graph_output = graph_output
  }

  return(output)
}


