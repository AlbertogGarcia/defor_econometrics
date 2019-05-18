### creating function to aggregate pixels into larger property, grid cell, county, etc. for analysis ###
library(tidyverse)


#county_fnc <- function(df, S){

### assign each pixel to a location on X-Y grid ###
# separate each unique pixel
unique <- unique(defor_df$idx)

idtreat <- subset(defor_df, select = c(idx, treat))
idtreat <- idtreat %>% distinct(idx, treat, .keep_all = TRUE)
idtreat$treat <- as.numeric(idtreat$treat)
idtreat <- idtreat[idtreat$treat == 1 , ]
iduntreat <- idtreat[idtreat$treat == 0 , ]
idtreat <- subset(idtreat, select = -c(treat))
iduntreat <- subset(iduntreat, select = -c(treat))
treatsamp <- sample_n(idtreat, nrow(idtreat))
untreatsamp <- sample_n(iduntreat, nrow(iduntreat))


# number of rows and columns for NxN grid

split(d, ceiling(seq_along(d)/20))

# NxN grid with random arrangement
grid <- matrix(rand, nrow)





### should have grid showing outcomes for each given year ###




### group into counties based on how many kilometers each should be ###

pixelper <- S/30

### calculate average deforestation rate for each county ###

#}
### end function ###






### use funvvtion to re-run analysis and cluster SEs at county level