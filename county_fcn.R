### creating function to aggregate pixels into larger property, grid cell, county, etc. for analysis ###
library(tidyverse)
library(gridExtra)

S=10
#county_fnc <- function(df, S){
y = S*1000
x <- y/30

### assign each pixel to a location on X-Y grid ###


# separate each unique pixel
unique <- unique(defor_df$idx)

idtreat <- subset(defor_df, select = c(idx, treat))
idtreat <- idtreat %>% distinct(idx, treat, .keep_all = TRUE)
idtreat$treat <- as.numeric(idtreat$treat)
idtreated <- idtreat[idtreat$treat == 0 , ]
iduntreated <- idtreat[idtreat$treat < 1 , ]
idtreated <- subset(idtreated, select = -c(treat))
iduntreated <- subset(iduntreated, select = -c(treat))
treatsamp <- sample_n(idtreated, nrow(idtreated))
untreatsamp <- sample_n(iduntreated, nrow(iduntreated))


# creating random groups of treated/untreated units

treatgroups <- split(treatsamp$idx, ceiling(seq_along(treatsamp$idx)/x))
untreatgroups <- split(untreatsamp$idx, ceiling(seq_along(untreatsamp$idx)/x))


#turning randomly assigned group vectors into county matrices
gridlist <- list()
for (i in 1:length(treatgroups)){
  treatmat <- matrix(unlist(treatgroups[i]), nrow = as.numeric(ceiling(sqrt(lengths(treatgroups[i])))), byrow = TRUE )
  nam <- paste("A", i, sep = "")
  assign(nam, treatmat)


  untreatmat <- matrix(unlist(untreatgroups[i]), nrow = as.numeric(ceiling(sqrt(lengths(untreatgroups[i])))), byrow = TRUE )
  name <- paste("B", i, sep = "")
  assign(name, untreatmat)
  
  gridlist <- c(gridlist, nam, name)
  
  
}
gridlist <- sample(gridlist)


### should have grid showing outcomes for each given year ###



### calculate average deforestation rate for each county ###

bob <- ceiling(sqrt(lengths(treatgroups[1])))


### reorient original dataframe

#}
### end function ###






### use funvvtion to re-run analysis and cluster SEs at county level