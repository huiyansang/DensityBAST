###################################################
#################### Libraries ####################
###################################################
rm(list=ls())
library(BASTION)
library(DensityBAST)
library(igraph)
library(dplyr)
library(sf)
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)
library(spatstat)
library(furrr)
library(ComplexDomain)
#####################
coords_mph=read.csv('coords_mph.csv')[,2:3]/1000
sf_coords=coords_mph %>%st_as_sf(coords = c("V1","V2"), crs = NA)
sf_bnd <- sf_coords%>%
  st_union() %>%
  st_convex_hull()%>%
  st_buffer(0.2)

training_coords=as.matrix(coords_mph);colnames(training_coords)=c('X1','X2')
plot(sf_bnd)
points(training_coords)

###
hyperpars = list(
    N_LEARNERS = 20,
    MESH_REF = 500,
    MCMC = 2000,
    BURNIN = 1000,
    parallel_cores = 8
  )

densityBASTfit = DensityBASTFit(sf_bnd,training_coords,hyperpars)


test_df=genVMesh(sf_bnd,n_ref=10000,type='regular',graph=FALSE)%>%st_as_sf()

test_df = test_df %>%
  mutate(est_intens = test_df%>%st_centroid()%>%st_coordinates()%>% densityBASTfit$bastion_estFunc_postmean())


ggplot() +
  geom_sf(data = test_df, aes(fill = log(est_intens)),color=NA,linewidth=0) +
  geom_point(data = data.frame(training_coords),
           aes(x = X1, y = X2), color = "#CC0000",
           shape = 8,
           size = 0.4,
           alpha = 0.25) +
  scale_fill_viridis_c()
