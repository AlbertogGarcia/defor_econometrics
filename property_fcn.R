# function to aggregate pixels to property level, adding in random perterbations at property level
library(sf)
library(rlist)

psize <- 1000
cellsize = 20



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



#### need to assign units location in treated or untreatedd gridcells

# randomly sampling locations in each type of grid
treat_locs <- st_sample(treat_grids, nobs)
untreat_locs <- st_sample(untreat_grids, nobs)#, coords = c("", ""))

#determine which county a point lies in
treat_counties <- st_within(treat_locs, treat_grids, sparse = FALSE, prepared = TRUE)*1
untreat_counties <- st_within(untreat_locs, untreat_grids, sparse = FALSE, prepared = TRUE)*1
t_whichcounty <- max.col(treat_counties)
u_whichcounty <- max.col(untreat_counties)

#determine which property a point lies in


plot(landscape)
plot(treat_locs, pch = 20, col = "red", add = TRUE)
plot(untreat_locs, pch = 20, col = "green", add = TRUE)
plot(overgrid, axes = TRUE, add = TRUE)
plot(v, axes = TRUE, add = TRUE)

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

### aggregate df to property level
suppressWarnings(
  propertylevel_df <-  aggregate(df_match, by = list(df_match$property, df_match$treat, df_match$year), FUN = mean, drop = TRUE)[c("property", "county", "treat", "year","defor")]
)

suppressWarnings(
  countylevel_df <-  aggregate(df_match, by = list(df_match$county, df_match$treat, df_match$year), FUN = mean, drop = TRUE)[c("county", "treat", "year","defor")]
)
countylevel_df <- unite(countylevel_df, county, treat, county, sep = "_", remove = FALSE)



# return new grid defor rate dataframe

# }

