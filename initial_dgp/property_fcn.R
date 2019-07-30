# function to aggregate pixels to property level, adding in random perterbations at property level
library(sf)
library(rlist)
library(tidyverse)

psize <- 1000
cellsize = 15

#starting function
property_fcn <- function(defor_df, psize, cellsize){

# turn dropped outcome observations to 1
defor_df$defor[is.na(defor_df$defor)] <- 1
  
nobs<- length(unique(defor_df$idx))
  
rootn <- ceiling(sqrt(nobs*2))
landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))

overgrid <- st_make_grid(landscape, cellsize, square = TRUE)

#trim grid to landscape
overgrid <- st_intersection(overgrid, landscape)

#now we want to randomly determine treated and untreated gridcells
treat_grids <- list.sample(overgrid, round(length(overgrid)/2))

# calculate untreated gridcells
untreat_grids <- overgrid[!overgrid %in% treat_grids]

# generate voronoi points
vorpts <- st_sample(landscape, psize)

v <- vorpts %>%  # consider the sampled points
  st_geometry() %>% #  as geometry only 
  st_union() %>% # unite them 
  st_voronoi() %>% # perform the voronoi tessellation
  st_collection_extract(type = "POLYGON") %>% # select the polygons
  st_intersection(overgrid)  # limit to within county boundaries

#getting areas for properies and counties
pareas <- data.frame(matrix(unlist(st_area(v))))
pareas <- tibble::rownames_to_column(pareas)
names(pareas)[1] <- paste("property")
names(pareas)[2] <- paste("parea")

t_careas <- data.frame(matrix(unlist(st_area(treat_grids))))
u_careas <- data.frame(matrix(unlist(st_area(untreat_grids))))
t_careas <- tibble::rownames_to_column(t_careas)
u_careas <- tibble::rownames_to_column(u_careas)
names(t_careas)[1] <- paste("county")
names(u_careas)[1] <- paste("county")
t_careas$treat <- rep(1, nrow(t_careas))
u_careas$treat <- rep(0, nrow(u_careas))
t_careas <- unite(t_careas, countyid, treat, county, sep = "_", remove = TRUE)
u_careas <- unite(u_careas, countyid, treat, county, sep = "_", remove = TRUE)
names(t_careas)[2] <- paste("carea")
names(u_careas)[2] <- paste("carea")
careas <- rbind(t_careas, u_careas)
# #### need to assign units location in treated or untreatedd gridcells

# randomly sampling locations in each type of grid
treat_locs <- st_sample(treat_grids, nobs)
untreat_locs <- st_sample(untreat_grids, nobs)#, coords = c("", ""))

#determine which county a point lies in
treat_counties <- st_within(treat_locs, treat_grids, sparse = FALSE, prepared = TRUE)*1
untreat_counties <- st_within(untreat_locs, untreat_grids, sparse = FALSE, prepared = TRUE)*1
t_whichcounty <- max.col(treat_counties)
u_whichcounty <- max.col(untreat_counties)


#### now we need to assign the random locations to the pixels in defor_df
# separating treated and untreated units
defor_df$treat <- as.numeric(defor_df$treat)
df_tr <- subset(defor_df, treat ==1)
df_un <- subset(defor_df, treat ==0)

# assigning each unique unit a location
unique_un <- unique(df_un$idx)
unique_tr <- unique(df_tr$idx)
df_un <- st_sf(unique_un, u_whichcounty, geometry = untreat_locs)
df_tr <- st_sf(unique_tr, t_whichcounty, geometry = treat_locs)

names(df_un)[2] <- paste("county")
names(df_tr)[2] <- paste("county")
colnames(df_tr)[1] <- "idx"
colnames(df_un)[1] <- "idx"
df_3 <- rbind(df_tr, df_un)

df_match <-  merge(df_3, defor_df, by = "idx")

whichprop <- st_within(df_match, v, sparse = FALSE, prepared = TRUE)*1
whichprop <- max.col(whichprop)
df_match$property <- whichprop

#match back in county and property areas merge

### aggregate df to property level
suppressWarnings(
  propertylevel_df <-  aggregate(df_match, by = list(df_match$property, df_match$treat, df_match$year), FUN = mean, drop = TRUE)[c("property", "county", "treat", "year","defor")]
)

suppressWarnings(
  countylevel_df <-  aggregate(df_match, by = list(df_match$county, df_match$treat, df_match$year), FUN = mean, drop = TRUE)[c("county", "treat", "year","defor")]
)

propertylevel_df <- unite(propertylevel_df, countyid, treat, county, sep = "_", remove = FALSE)
propertylevel_df <-  merge(propertylevel_df, pareas, by = "property")
countylevel_df <- unite(countylevel_df, countyid, treat, county, sep = "_", remove = FALSE)
countylevel_df<-  merge(countylevel_df, careas, by = "countyid")

assign('countylevel_df',countylevel_df, envir=.GlobalEnv)
assign('propertylevel_df',propertylevel_df, envir=.GlobalEnv)

# return new grid defor rate dataframe

}

