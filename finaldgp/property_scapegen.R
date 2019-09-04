# function to generate generic landscape
# pixel value realizations can then be simulated afterward
library(sf)
library(rlist)
library(tidyverse)

property_scapegen <- function(nobs, cellsize, ppoints){
  
  
  rootn <- ceiling(sqrt(nobs))
  landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))
  
  overgrid <- st_make_grid(landscape, cellsize, square = TRUE)
  
  #trim grid to landscape
  overgrid <- st_intersection(overgrid, landscape)
  
  
  #generate pixels within landscape
  DT <- data.frame(
    pixels= seq(from = 1, to = rootn^2),
    longitude= (rep(seq(from = 1, to = rootn), rootn) - .33),
    latitude= (rep(1:rootn, each= rootn) - .33)
  ) 
  
  pixloc_df <- st_as_sf(DT, coords = c("longitude", "latitude"))
  
  #generate voronoi pts for 
  vorpts <- st_sample(landscape, ppoints)
  
  v <- vorpts %>%  # consider the sampled points
    st_geometry() %>% #  as geometry only 
    st_union() %>% # unite them 
    st_voronoi() %>% # perform the voronoi tessellation
    st_collection_extract(type = "POLYGON") %>% # select the polygons
    st_intersection(overgrid)  # limit to within grid boundaries
  
  #determine which pixels are in each grid 
  wgrid <- st_within(pixloc_df, overgrid, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$grid <- max.col(wgrid)
  
  wprop <- st_within(pixloc_df, v, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$property <- max.col(wprop)
  
  #determine treated vs. untreated grids
  treat_grids <- sample(1:length(overgrid), (length(overgrid)/2) )
  grid= seq(from = 1, to = length(overgrid))
  treatgrid <- data.frame(
    grid= seq(from = 1, to = length(overgrid)),
    treat = ifelse(grid %in% treat_grids, 1, 0)
  ) 
  
  # determine which pixels are treated vs. untreated
  pixloc_df <- merge(pixloc_df, treatgrid, by = "grid")
  
  #getting areas for properies and counties
  pareas <- data.frame(matrix(unlist(st_area(v))))
  pareas <- tibble::rownames_to_column(pareas)
  names(pareas)[1] <- paste("property")
  names(pareas)[2] <- paste("parea")
  gareas <- data.frame(matrix(unlist(st_area(overgrid))))
  gareas <- tibble::rownames_to_column(gareas)
  names(gareas)[1] <- paste("grid")
  names(gareas)[2] <- paste("garea")
  
  #merging back areas
  areas <- merge(pixloc_df, gareas, by = "grid")
  pixloc_df <- merge(areas, pareas, by = "property")
  
  outputs = list('pixloc_df' = pixloc_df)
  return(outputs)
  
}
