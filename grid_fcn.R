### the purpose of this function is to aggregate pixels to a county or property level
library(sp)
library(sf)
library(ggplot2)
library(rlist)
library(tidyverse)

# grid_fcn <- function(defor_df, celllsize) {

 

# first, generate landscape based on size
# there are nobs treated and untreated obs for nobs*2 total obs in each period
rootn <- ceiling(sqrt(nobs*2))
landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))

overgrid <- st_make_grid(landscape, cellsize = 11, square = TRUE)

#trim grid to landscape
overgrid <- st_intersection(overgrid, landscape)

#now we want to randomly determine treated and untreated gridcells
treat_grids <- list.sample(overgrid, round(length(overgrid)/2))

# calculate untreated gridcells
untreat_grids <- setdiff(overgrid, treat_grids)

# could plot treated cells over landscape
# plot(landscape, col = "#339900")
# plot(overgrid, axes = TRUE, add = TRUE)
# plot(treat_grids, col = "red", add = TRUE)

# need to assign units location in treated or untreatedd gridcells

st_sample(overgrid, 100)


plot(landscape, col = "#339900")
plot(overgrid, axes = TRUE, add = TRUE)
plot(st_sample(overgrid, 100), col = "red", add = TRUE)









