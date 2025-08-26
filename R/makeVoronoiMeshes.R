#' Generate list of random Voronoi meshes on a boundary
#'
#' @param boundary
#' @param N_LEARNERS
#' @param MESH_REF
#'
#' @return A list of length N_LEARNERS where each element is a collection of sf polygons
#' @export
#'
#' @examples
makeVoronoiMeshes = function(boundary, N_LEARNERS, MESH_REF, max_area_ratio = NULL) {

  custom_mesh = vector("list", N_LEARNERS)

  for(learner in 1:N_LEARNERS) {
    # if(data_mesh) {
    #   custom_mesh[[learner]] = sample(which(partition[[learner]] != learner),
    #                                   size = MESH_REF,
    #                                   replace = TRUE) %>%
    #     data[.] %>%
    #     do.call(c, .) %>%
    #     sf::st_voronoi() %>%
    #     sf::st_collection_extract() %>%
    #     sf::st_intersection(boundary) %>%
    #     sf::st_sf() %>%
    #     sf::st_cast("MULTIPOLYGON") %>%
    #     sf::st_cast("POLYGON")
    # } else {
      custom_mesh[[learner]] = sf::st_sample(boundary, size = MESH_REF) %>%
        do.call(c, .) %>%
        sf::st_voronoi() %>%
        sf::st_collection_extract() %>%
        sf::st_intersection(boundary) %>%
        sf::st_sf() %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_cast("POLYGON")
    # }
  }
  # Optionally, ensure that area of largest polygon is not more than
  # max_area_ratio*area of smallest polygon
  if(!is.null(max_area_ratio)) {
    for(learner in 1:N_LEARNERS) {
      areas = sf::st_area(custom_mesh[[learner]])
      area_ratio = max(areas)/min(areas)
      while(area_ratio > max_area_ratio) {
       # Until condition is satisfied, merge smallest polygon with the
       # smallest of its neighbors (alternative idea is to merge along largest border)
        smallest_polygon_idx = which.min(areas)
        # Neighbors
        merging_neighbor_idx = sf::st_relate(custom_mesh[[learner]][smallest_polygon_idx,],
                                              custom_mesh[[learner]][-smallest_polygon_idx,],
                                              pattern = "****1****")[[1]] %>%
                                              # Smallest of the neighbors
                                              `[`(., {`[`(areas, .) %>% which.min})
        # Get merged polygon
        merged_polygon = sf::st_union(custom_mesh[[learner]][smallest_polygon_idx,],
                                  custom_mesh[[learner]][merging_neighbor_idx,]) %>%
                          sf::st_cast("POLYGON") %>%
                          sf::st_sf()
        # Update custom_mesh[[learner]] by removing the two and adding the merged
        #custom_mesh[[learner]][smallest_polygon_idx,] = merged_polygon
        custom_mesh[[learner]] = custom_mesh[[learner]][-c(merging_neighbor_idx, smallest_polygon_idx),]
        custom_mesh[[learner]] = rbind(custom_mesh[[learner]], merged_polygon)

        areas = sf::st_area(custom_mesh[[learner]])
        area_ratio = max(areas)/min(areas)
      }
    }
  }

  return(custom_mesh)
}
