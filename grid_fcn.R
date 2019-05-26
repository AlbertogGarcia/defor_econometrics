### the purpose of this function is to aggregate pixels to a county or property level
library(sp)
library(sf)
library(ggplot2)
library(rlist)
library(tidyverse)

cellsize = 11

# grid_fcn <- function(defor_df, cellsize){

 

# first, generate landscape based on size
# there are nobs treated and untreated obs for nobs*2 total obs in each period
rootn <- ceiling(sqrt(nobs*2))
landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))

overgrid <- st_make_grid(landscape, cellsize, square = TRUE)

#trim grid to landscape
overgrid <- st_intersection(overgrid, landscape)

#now we want to randomly determine treated and untreated gridcells
treat_grids <- list.sample(overgrid, round(length(overgrid)/2))

# calculate untreated gridcells
untreat_grids <- overgrid[!overgrid %in% treat_grids]


#### need to assign units location in treated or untreatedd gridcells

# randomly sampling locations in each type of grid
treat_locs <- st_sample(treat_grids, nobs)
untreat_locs <- st_sample(untreat_grids, nobs)


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

#matching coordinates back to original dataframe

df_match <-  merge(df_3, defor_df, by = "idx")


### calculate average defor rate in each county based on year
# aggregate up to county in each year 
suppressWarnings(
  countylevel_df <-  aggregate(df_match, by = list(df_match$county, df_match$treat, df_match$year), FUN = mean, drop = TRUE)[c("county", "treat", "year","defor")]
)
countylevel_df <- unite(countylevel_df, county, treat, county, sep = "_", remove = FALSE)


# }

