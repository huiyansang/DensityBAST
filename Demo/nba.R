library(sf)
library(DensityBAST)
library(ggplot2)
load("sf_pp_curry.Rdata")

thinning_p = 0.8
N = nrow(coords)
# "training" set
point_included = purrr::rbernoulli(N, p = thinning_p)
training_coords = coords[point_included, ]
training_sf_coords = sf_coords[point_included, ]
# "testing" set
testing_coords = coords[!point_included, ]
testing_sf_coords = sf_coords[!point_included, ]


training_points = training_coords %>%
  st_multipoint() %>%
  st_sfc %>%
  st_cast(to = "POINT")
# densityBAST fit
densityBASTfit = DensityBASTFit(sf_bnd,training_coords)

# horseshoe example
#load('horseshoe.Rdata')
#densityBASTfit = DensityBASTFit(horseshoe.domain,horseshoe.pp)


## plot the estimated intensity
eval_tesselation = st_make_grid(sf_bnd, cellsize = c(1, 1))

grid_df = st_make_grid(sf_bnd, cellsize = c(0.25, 0.25)) %>%
  st_centroid %>%
  st_coordinates %>%
  data.frame()

plot_df = grid_df %>%
  rename(x = 1, y = 2) %>%
  mutate(est_intens = grid_df[,1:2] %>% as.matrix %>% (densityBASTfit$bastion_estFunc_postmean))

estimated_mae_yin = function(p,
                             tesselation, # uniform squares in Yin et al.
                             holdout_points, #sf coords
                             est_intensity_func) # function of coords
{

  centroids = tesselation %>%
    st_centroid %>%
    st_coordinates %>%
    `colnames<-`(c('x','y'))
  areas = tesselation %>% st_area
  est_intens = centroids %>% est_intensity_func
  # assumes locally homogeneous estimated intensity
  # over each polygon in the tesselation
  expected_num_points = ((1-p)/p)*(areas * est_intens)
  observed_num_points = holdout_points %>%
    st_contains(tesselation, .) %>%
    sapply(., length)

  abs_error = abs(expected_num_points - observed_num_points)
  mae = mean(abs_error)
  return(mae)
}


nba_example_mae = estimated_mae_yin(thinning_p, eval_tesselation,
                                    testing_sf_coords, densityBASTfit$bastion_estFunc_postmean)

nba_plot_court = ggplot() +
  geom_raster(data = plot_df, aes(x = x, y = y, fill = est_intens)) +
  geom_point(data = data.frame(training_coords),
             aes(x = X1, y = X2), shape = 8, size = 1.5, alpha = 0.3) +
  scale_fill_gradient(low = "#FFFFFF", high = "#003C71", trans = "log1p") +
  #geom_sf(data = court_sf, fill = NA) +
  coord_fixed() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") + ylab("") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())  +
  theme_void() +
  theme(legend.position = "none")

nba_plot_court
