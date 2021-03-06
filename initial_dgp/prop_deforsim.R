# adding in property level perturbations in order to see which aggregation method is best

#packages
library(dplyr)
library(purrr)
library(tidyverse)
library(reshape2)
library(data.table)
library(sf)
library(rlist)

#define parameters
years <- 3               # number of years in each of two periods
nobs <- 10000           # number of observations in each of two groups
psize <- 1000
cellsize = 25

#starting function
prop_deforsim <- function(nobs, years, psize, cellsize){

d_outside <- 0.30        # true deforestation rate outside treated area in period 2 
diff <- 0.40             # pre-treatment difference between treatment and control
trend <- -0.10           # trend in deforestation rate across periods
att <- -0.16             # average treatement effect on the treated


#deforestation rates
m_00 <- d_outside                        # average deforestation rate outside treated area in p1
m_01 <- d_outside + trend                # average deforestation rate outside treated area in p2
m_10 <- d_outside + diff                 # average deforestation rate inside treated area in p1
m_11 <- d_outside + diff + trend + att   # average deforestation rate inside treated area in p2


# introducing three types of randomness
#     i: Individual-level
#     y: year-level changes
#     e: common structure across all observations

std_i <- 0.1
std_y <- 0.01
std_e <- 0.01


e_err <- matrix(rnorm(nobs*2 * years*2 , 0 , std_e), ncol = years*2)
i_err <- replicate(years*2, rnorm(nobs*2 , 0 , std_i))
y_err <- t(replicate(nobs*2, rnorm(years*2 , 0 , std_y)))

df <- data.frame( i_err + y_err + e_err)


#Add in average deforestation rates for each group
df[1:nobs,1:years] <- df[1:nobs,1:years] + m_00
df[1:nobs,(years+1):ncol(df)] <- df[1:nobs,(years+1):ncol(df)] + m_01
df[(nobs+1):nrow(df),1:years] <- df[(nobs+1):nrow(df),1:years] + m_10
df[(nobs+1):nrow(df),(years+1):ncol(df)] <- df[(nobs+1):nrow(df),(years+1):ncol(df)] + m_11


#define treatment
df$treat <- rep(0, nrow(df))
df$treat[(nobs+1):nrow(df)] <- df$treat[(nobs+1):nrow(df)]+1



##################################################################
####### Adding property level perturbations ######################
# up to this point we have each obs and value for four years


rootn <- ceiling(sqrt(nobs*2))
landscape = st_sfc(st_polygon(list(cbind(c(0,rootn,rootn,0,0),c(0,0,rootn,rootn,0)))))


#create grid and trim to landscape
overgrid <- st_intersection(st_make_grid(landscape, cellsize, square = TRUE), landscape)

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

# generating perturbations for each property
std_p <- 0.5
p_err <- data.frame(rnorm(length(v) , 0 , std_p))
p_err <- tibble::rownames_to_column(p_err)
names(p_err)[1] <- paste("property")

#### need to assign units location in treated or untreatedd gridcells

# randomly sampling locations in each type of grid
treat_locs <- st_sample(treat_grids, nobs)
untreat_locs <- st_sample(untreat_grids, nobs)#, coords = c("", ""))

#determine which county a point lies in
t_whichcounty <- max.col(st_within(treat_locs, treat_grids, sparse = FALSE, prepared = TRUE)*1)
u_whichcounty <- max.col(st_within(untreat_locs, untreat_grids, sparse = FALSE, prepared = TRUE)*1)


#### now we need to assign the random locations to the pixels in defor_df

df$treat <- as.numeric(df$treat)

# separating treated and untreated units and assigning each unique unit a location
df_un <- st_sf(subset(df, treat ==0), u_whichcounty, geometry = untreat_locs)
df_tr <- st_sf(subset(df, treat ==1), t_whichcounty, geometry = treat_locs)


names(df_un)[(years*2+2)] <- paste("county")
names(df_tr)[(years*2+2)] <- paste("county")
# colnames(df_tr)[1] <- "idx"
# colnames(df_un)[1] <- "idx"
df <- rbind(df_tr, df_un)

whichprop <- st_within(df, v, sparse = FALSE, prepared = TRUE)*1
whichprop <- max.col(whichprop)
df$property <- whichprop
names(p_err)[2]<-"p_err"


# df <-  merge(df, p_err, by = "property")
df <- tibble::rownames_to_column(merge(df, p_err, by = "property"))

df <- unite(df, idx, treat, rowname, sep = "_", remove = FALSE)


#adding property level errors
rownames(df) <- df$idx
df <- subset(df, select = -c(rowname, idx))




df_comp <- df[,c(2:(years*2+1))] + df$p_err


df_comp <- subset(df_comp, select = -c(geometry))


df <- subset(df, select = c(property, county, treat))
df <- tibble::rownames_to_column(df)



###################################################################




### simulating deforestation ###
defor_draw <- matrix(runif(nobs*2 * years*2 , 0 , 1), ncol = 4)
defor_df <- data.frame( df_comp > defor_draw)*1
not_defor <- rowSums(defor_df)<1 *1
defor_year <- max.col(defor_df, ties.method = "first")       #creating defor_year variable


defor_df <- cbind(defor_df, defor_year)
defor_df <- transform(defor_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))


defor_df <- tibble::rownames_to_column(defor_df)


#names(defor_df)[2:(years*2+1)] <- c(1:(years*2))

#defor_df <-  merge(df, defor_df, by = "rowname")
#defor_df <- subset(defor_df, select = -c(county))






defor_df <- melt(defor_df, id.vars = c('rowname', 'defor_year'))


setnames(defor_df, old=c("variable","value"), new=c("year", "defor"))

#replace defor with NA for any year > defor_year
defor_df$year <- as.numeric(defor_df$year)
defor_df$indic <- (defor_df$year - defor_df$defor_year)
#defor_df$defor <- ifelse(defor_df$indic > 0 , NA, defor_df$defor)
defor_df$defor <- ifelse(defor_df$indic > 0 , 1, defor_df$defor)
defor_df <- subset(defor_df, select = -c(indic))


defor_df$post <- (defor_df$year > years)*1


defor_df <-  merge(df, defor_df, by = "rowname")
rownames(df) <- df$rowname

defor_df <- subset(defor_df, select = -c(rowname))

suppressWarnings(
  proppert_df <-  aggregate(defor_df, by = list(defor_df$property, defor_df$treat, defor_df$year), FUN = mean, drop = TRUE)[c("property", "post", "county", "treat", "year","defor")]
)

proppert_df <- unite(proppert_df, countyid, treat, county, sep = "_", remove = FALSE)

suppressWarnings(
  countypert_df <-  aggregate(defor_df, by = list(defor_df$county, defor_df$treat, defor_df$year), FUN = mean, drop = TRUE)[c("county", "treat", "year","defor", "post")]
)
countypert_df <- unite(countypert_df, countyid, treat, county, sep = "_", remove = FALSE)
countypert_df <- subset(countypert_df, select = -c(county))


assign('proppert_df',proppert_df, envir=.GlobalEnv)

assign('countypert_df',countypert_df, envir=.GlobalEnv)

# end function
}
