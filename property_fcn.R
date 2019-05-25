# function to aggregate pixels to property level, adding in random perterbations at property level
library(sf)
library(MASS)
library(deldir)
library(muHVT)


# proposed layout of function:
# 1) create polygon with size of landscape
# 2) generate random points in polygon and perform voneroi tesselation to determine county boundaries
# 3) randomly assign counties to treatment or control status
# 4) randomly assign units to location on landscape depending on treatment status
# 5) calculate average defor rate within each county
# 6) 

# #now within the landscape, we generate 10 points from which to generate our voronoi tesselations
# points = 10
# 
# x <- sample(1:rootn,points)
# y <- sample(1:rootn,points)
# 
# pol <- data.frame(x, y)
# ggplot(pol,aes(x,y)) +
#   stat_voronoi(geom="path") +
#   geom_point()
psize <- 5000
rootn <- ceiling(sqrt(nobs*2))
landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))

overgrid <- st_make_grid(landscape, cellsize = 25, square = TRUE)

#trim grid to landscape
overgrid <- st_intersection(overgrid, landscape)

#now we want to randomly determine treated and untreated gridcells
treat_grids <- list.sample(overgrid, round(length(overgrid)/2))

# calculate untreated gridcells
untreat_grids <- overgrid[!overgrid %in% treat_grids]

# generate voronoi points
vor_pts <- st_sample(landscape, psize)
v <- vor_pts %>%  # consider the master points
  st_geometry() %>% # ... as geometry only (= throw away the data items)
  st_union() %>% # unite them ...
  st_voronoi() %>% # ... and perform the voronoi tessellation
  st_collection_extract(type = "POLYGON") %>% # select the polygons
  st_intersection(landscape)  # limit to Prague city boundaries


#### need to assign units location in treated or untreatedd gridcells

# randomly sampling locations in each type of grid
treat_locs <- st_sample(treat_grids, nobs)
untreat_locs <- st_sample(untreat_grids, nobs)#, coords = c("", ""))

plot(landscape)
plot(treat_locs, pch = 20, col = "brown", add = TRUE)
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
df_tr <- st_sf(unique_un, geometry = untreat_locs)
df_un <- st_sf(unique_tr, geometry = treat_locs)


# matching coordinates back to original dataframe
colnames(df_tr)[1] <- "idx"
colnames(df_un)[1] <- "idx"
df_3 <- rbind(df_tr, df_un)

df_match <-  merge(df_3, defor_df, by = "idx")



### calculate average defor rate in each county based on year



# return new grid defor rate dataframe

# }

