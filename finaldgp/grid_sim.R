library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
source('gridscapegen.R')

#begin function
grid_sim <- function(n, nobs, years, b0, b1, b2, b3, std_a = 0.1, std_v = 0.25, cellsize){

gridscape <- gridscapegen(nobs, cellsize)
pixloc_df = gridscape$pixloc_df
gridcoords = gridscape$gridcoords
N_treat = gridscape$N_treat    

pixloc <- pixloc_df[order(pixloc_df$pixels),]
  
coeffmatrix <- matrix(nrow = n, ncol = 1)

  for(i in 1:n){
    
    panels <- fabricate(
      pixels = add_level(N = nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
      year = add_level(N = (years*2), nest = FALSE),
      obs = cross_levels(
        by = join(pixels, year),
        post = ifelse(year > years, 1, 0),
        v_it = rnorm(N, 0, std_v),
        ystar = b0 + b1*treat + b2*post + b3*treat*post + a_i + v_it,
        y = (ystar > 0)*1
      )
    )
    #need to determine which year deforestation occurred
    year_df <- subset(panels, select = c(pixels, year, y))
    #year_df <- melt(year_df, id.vars = c("pixels", "y_it"), value.name = "year")
    year_df <- dcast(year_df, pixels ~ year , value.var = "y")
    rownames(year_df) <- year_df$pixels
    year_df <- subset(year_df, select = -c(pixels))
    
    #creating variable for the year a pixel is deforested
    not_defor <- rowSums(year_df)<1 *1
    defor_year <- max.col(year_df, ties.method = "first") 
    defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
    defor_df <- tibble::rownames_to_column(defor_df)
    names(defor_df)[1] <- paste("pixels")
    defor_df <- subset(defor_df, select = c(pixels, defor_year))
    panels <- merge(defor_df, panels, by = "pixels")
    
    
    # creating three outcome variables for each possible situation
    ### y: allows the outcome to switch between 0 and 1 across years
    ### y_it: outcome is dropped in years after pixel is first deforested
    ### defor: outcome is set to 1 in each year after the pixel is deforested
    panels$year <- as.numeric(panels$year)
    panels$indic <- (panels$year - panels$defor_year)
    panels$defor <- ifelse(panels$indic > 0 , 1, panels$y)
    panels <- subset(panels, select = -c(indic))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    
    panels <- merge(pixloc, panels, by = c("pixels", "treat"))
    
    
    # aggregate up to county in each year 
    suppressWarnings(
      gridlevel_df <-  aggregate(panels, by = list(panels$grid, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("grid", "treat", "post", "year","defor")]
    )
    
    return(gridlevel_df)
    
    
  }
  
  
  



}
