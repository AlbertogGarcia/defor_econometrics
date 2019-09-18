# function to generate generic landscape
# pixel value realizations can then be simulated afterward
library(sf)
library(rlist)
library(tidyverse)

county_scapegen <- function(nobs, cellsize, ppoints, cpoints){
  
  
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
  vorpts_prop <- st_sample(landscape, ppoints)
  vorpts_county <- st_sample(landscape, cpoints)
  
  v_county <- vorpts_county %>%  # consider the sampled points
    st_geometry() %>% #  as geometry only 
    st_union() %>% # unite them 
    st_voronoi() %>% # perform the voronoi tessellation
    st_collection_extract(type = "POLYGON")
  
  v_property <- vorpts_prop %>%  # consider the sampled points
    st_geometry() %>% #  as geometry only 
    st_union() %>% # unite them 
    st_voronoi() %>% # perform the voronoi tessellation
    st_collection_extract(type = "POLYGON") %>% # select the polygons
    st_intersection(v_county)  # limit to within county boundaries
  
  #determine which pixels are in each grid 
  wgrid <- st_within(pixloc_df, overgrid, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$grid <- max.col(wgrid)
  
  wprop <- st_within(pixloc_df, v_property, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$property <- max.col(wprop)
  
  wcounty <- st_within(pixloc_df, v_county, sparse = FALSE, prepared = TRUE)*1
  pixloc_df$county <- max.col(wcounty)
  
  #determine treated vs. untreated counties
  treat_counties <- sample(1:length(v_county), (length(v_county)/2) )
  county= seq(from = 1, to = length(v_county))
  
  treatcounty <- data.frame(
    county= seq(from = 1, to = length(v_county)),
    treat = ifelse(county %in% treat_counties, 1, 0)
  ) 
  
  # determine which pixels are treated vs. untreated
  pixloc_df <- merge(pixloc_df, treatcounty, by = "county")
  
  #getting areas for properies and counties
  pareas <- data.frame(matrix(unlist(st_area(v_property))))
  pareas <- tibble::rownames_to_column(pareas)
  names(pareas)[1] <- paste("property")
  names(pareas)[2] <- paste("parea")
  careas <- data.frame(matrix(unlist(st_area(v_county))))
  careas <- tibble::rownames_to_column(careas)
  names(careas)[1] <- paste("county")
  names(careas)[2] <- paste("carea")
  gareas <- data.frame(matrix(unlist(st_area(overgrid))))
  gareas <- tibble::rownames_to_column(gareas)
  names(gareas)[1] <- paste("grid")
  names(gareas)[2] <- paste("garea")
  
  #merging back areas
  areas <- merge(pixloc_df, gareas, by = "grid")
  c_areas <- merge(areas, careas, by = "county")
  pixloc_df <- merge(c_areas, pareas, by = "property")
  
  outputs = list('pixloc_df' = pixloc_df)
  return(outputs)
  
}