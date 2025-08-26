
#' Generate learners for fitting densityBAST
#'
#' @param meshes List of meshes, as returned by makeVoronoiMeshes or makeDelauneyMeshes
#' @param data An object of class sfc_POINT; a list of the coordinates of observed data locations
#'
#' @return
#' @export
#'
#' @examples
makeLearners = function(meshes, data) {

  N = length(data)
  N_LEARNERS = length(meshes)

  #meshes = vector("list", N_LEARNERS)
  full_graphs = vector("list", N_LEARNERS)
  polygon_areas = vector("list", N_LEARNERS)
  point_in_polygon = vector("list", N_LEARNERS)
  points_per_polygon = vector("list", N_LEARNERS)
  data_polygon_id = vector("list", N_LEARNERS)
  n_polygons = vector("list", N_LEARNERS)

  pb <- progress::progress_bar$new(
    format = "  Building Learners [:bar] :percent eta: :eta",
    total = N_LEARNERS, clear = FALSE, width= 60)

  for(learner in 1:N_LEARNERS) {
    pb$tick()
    full_graph = meshes[[learner]] %>%
      sf::st_relate(pattern = "****1****") %>%# Link each triangle with it's neighbors that share an edge
      igraph::graph.adjlist() %>%                 # Use these links to create a graph
      igraph::as.undirected() %>%                 # Make the graph undirected
      igraph::simplify()                  # Remove self edges / loops
    igraph::E(full_graph)$weight = runif(igraph::gsize(full_graph), 0, 1) # Put weights on the graph
    full_graphs[[learner]] = full_graph
    # Get the area of each polygon
    polygon_areas[[learner]] = sf::st_area(meshes[[learner]])
    # Get which data points lie in each polygon
    point_in_polygon[[learner]] = sf::st_contains(meshes[[learner]], data)
    # Get the number of points in each polygon
    points_per_polygon[[learner]] = sapply(point_in_polygon[[learner]], length)
    # For each data point, the polygon id the data point lies in
    data_polygon_id[[learner]] = get_polygon_id(point_in_polygon[[learner]])
    # Number of polygons in each learner
    n_polygons[[learner]] = length(polygon_areas[[learner]])
  }
  return(
    list(
      data = data,
      meshes = meshes,
      full_graphs = full_graphs,
      node_measures = polygon_areas,
      node_contains = point_in_polygon,
      node_counts = points_per_polygon,
      data_membership = data_polygon_id,
      graph_sizes = n_polygons
    )
  )

}
