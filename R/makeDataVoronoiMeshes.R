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
makeDataVoronoiMeshes = function(boundary, data_matrix,N_LEARNERS, MESH_REF, max_area_ratio = NULL) {

  custom_mesh = vector("list", N_LEARNERS)

  for(learner in 1:N_LEARNERS) {
   set.seed(learner+1234)
   data_matrix[sample(1:nrow(data_matrix), MESH_REF, replace = TRUE),] %>%
      st_multipoint() %>%
      st_sfc %>%
      st_cast(to = "POINT") -> coords_ref

    coords_ref = c(coords_ref,st_sample(sf_bnd, size = 50,'regular'))
    voronoi <- st_voronoi(st_union(coords_ref),sf_bnd)
    sf_mesh = st_intersection(st_cast(voronoi), sf_bnd)
    custom_mesh[[learner]] = sf_mesh
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
