library(DeclareDesign)
library(reshape2)  
deforestation_DGP <- function(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v){
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2)^.5)
  
  panels <- fabricate(
    pixels = add_level(N = nobs, a_i = rnorm(N, 0, std_a), treat = rbinom(N, 1, 0.5)),
    year = add_level(N = (years*2), nest = FALSE),
    obs = cross_levels(
      by = join(pixels, year),
      post = ifelse(year > years, 1, 0),
      v_it = rnorm(N, 0, std_v),
      ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
      y = (ystar > 0)*1
    )
  )
  
  #need to determine which year deforestation occurred
  year_df <- subset(panels, select = c(pixels, year, y))
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
  panels$y_it <- ifelse(panels$indic > 0 , NA, panels$y)
  panels$defor <- ifelse(panels$indic > 0 , 1, panels$y_it)
  panels <- subset(panels, select = -c(indic))
  
  outputs = list("panels" = panels, "ATT" = ATT)
  # assign('panels',panels, envir=.GlobalEnv)
  # assign('ATT',ATT, envir=.GlobalEnv)
  return(outputs)
}


