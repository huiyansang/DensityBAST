
#' Generate list of random Delauney Triangulation meshes on a boundary
#'
#' @param boundary
#' @param N_LEARNERS
#' @param MESH_REF
#' @param MAX_ATTEMPTS
#' @param eps
#'
#' @return
#' @export
#'
#' @examples
makeDelauneyMeshes = function(boundary,
                                N_LEARNERS,
                                MESH_REF,
                                MAX_ATTEMPTS = 100,
                                eps = 0.0000001) {

  custom_mesh = vector("list", N_LEARNERS)

  for(learner in 1:N_LEARNERS) {
    attempts = 0
    repeat {
      custom_mesh[[learner]] = boundary %>% genMesh(n_ref = MESH_REF)
      polygon_areas = sf::st_area(custom_mesh[[learner]])
      attempts = attempts + 1

      if(all(polygon_areas >= eps)) {
        break
      }
      if(attempts > MAX_ATTEMPTS) {
        warning(paste0("Weak learner ", learner, " triangulation: Triangle areas too small. Consider changing MESH_REF\n"))
        break
      }
    }
  }
  return(custom_mesh)
}
