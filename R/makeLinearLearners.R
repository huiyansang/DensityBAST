
getLintessNeighbors = function(lintess) {
  split_seg_nums = lintess$df$seg %>% table() %>% `>=`(., 2) %>% which()
  split_seg_ids = which(lintess$df$seg %in% split_seg_nums)

  split_segs = lintess$df[split_seg_ids,]

  n_factors = length(levels(lintess$df$tile))
  output_adj_mat = matrix(0, nrow = n_factors, ncol = n_factors)

  for(i in 1:(nrow(split_segs) - 1)) {
    seg0 = split_segs$seg[i]
    seg1 = split_segs$seg[i+1]
    tile0 = split_segs$tile[i]
    tile1 = split_segs$tile[i+1]
    t1 = split_segs$t1[i]
    t0_next = split_segs$t0[i+1]
    # 1 if true, 0 if false
    output_adj_mat[tile0, tile1] = as.numeric(
      (seg0 == seg1)*(tile0 != tile1)*(t1 == t0_next) || output_adj_mat[tile0, tile1]
    )
  }
  return(output_adj_mat)
}

#' Title
#'
#' @param linenetwork_lpp A spatstat.linnet::lpp class containing the observed points and line geometry
#' @param N_LEARNERS Number of weak learners
#' @param MESH_REF Number of mesh reference nodes per weak learner
#' @param usepoisson Boolean, whether to simulate reference nodes from a poisson point process instead of uniformly
#'
#' @return
#' @export
#'
#' @examples
makeLinearLearners = function(data_matrix,
                              N_LEARNERS,
                              MESH_REF,
                              usepoisson = TRUE) {
  # Original input containing both the domain information (a linnet)
  # and the observed points (a psp)
  data = data_matrix
  N = nrow(data$data)
  networkLength = volume(data_matrix$domain)

  meshes = vector("list", N_LEARNERS)
  full_graphs = vector("list", N_LEARNERS)
  tile_lengths = vector("list", N_LEARNERS)
  point_in_tile = vector("list", N_LEARNERS)
  points_per_tile = vector("list", N_LEARNERS)
  data_tile_id = vector("list", N_LEARNERS)
  n_tiles = vector("list", N_LEARNERS)

  pb <- progress::progress_bar$new(
    format = "  Building Learners [:bar] :percent eta: :eta",
    total = N_LEARNERS, clear = FALSE, width= 60)


  # Build meshes
  for(learner in 1:N_LEARNERS) {
    pb$tick()
    if(usepoisson) {
      mesh_lpp = spatstat.linnet::rpoislpp(MESH_REF/(networkLength), as.linnet(data), drop = TRUE, nsim = 1)
    } else {
      # Uniformly sample reference points on the linnet
      mesh_lpp = spatstat.linnet::runiflpp(MESH_REF, as.linnet(data), drop = TRUE, nsim = 1)
    }
    # Perform voronoi tesselation based on sampled points
    tesselation = spatstat.linnet::lineardirichlet(mesh_lpp)
    meshes[[learner]] = tesselation
    # Get lengths of each tile (areas of each polygon)
    tile_lengths[[learner]] = spatstat.linnet::tile.lengths(tesselation)
    # Build neighbor graph
    full_graph = getLintessNeighbors(tesselation) %>% # Link each tile with it's neighbors that share an edge
      igraph::graph_from_adjacency_matrix() %>%                 # Use these links to create a graph
      igraph::as_undirected() %>%                 # Make the graph undirected
      igraph::simplify()                  # Remove self edges / loops
    igraph::E(full_graph)$weight = runif(igraph::gsize(full_graph), 0, 1)
    full_graphs[[learner]] = full_graph
    # Verify the correct number of nodes in the graph
    n_tiles[[learner]] = spatstat.linnet::nobjects.lintess(tesselation)
    # Bin observed points into tesselation
    # For each data point, the tile id the data point lies in
    data_tile_id[[learner]] = spatstat.linnet::as.linfun.lintess(tesselation)(data) %>% as.numeric()
    # Get which data points lie in each tile
    point_in_tile[[learner]] = lapply(1:n_tiles[[learner]], function(x){which(data_tile_id[[learner]] == x)})
    # Get the number of points in each tile
    points_per_tile[[learner]] = sapply(point_in_tile[[learner]], length)
  }

  return(
    list(
      data = data,
      meshes = meshes,
      full_graphs = full_graphs,
      node_measures = tile_lengths,
      node_contains = point_in_tile,
      node_counts = points_per_tile,
      data_membership = data_tile_id,
      graph_sizes = n_tiles
    )
  )

}
