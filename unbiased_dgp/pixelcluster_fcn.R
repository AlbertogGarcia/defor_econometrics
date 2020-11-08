#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level


library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tidyverse)
library(tictoc)
source('full_landscape.R')

#begin function
pixelcluster_fcn <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, std_p = 0.0, cellsize, ppoints, cpoints){
  
  countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  pixloc <- pixloc_df
  
  covermat <- matrix(nrow = n, ncol = 1)
  coeffmatrix <- matrix(nrow = n, ncol = 1)
  
  for(i in 1:n){
    tic("loop")
    Nobs <- length(pixloc$treat)  
    panels <- fabricate(
      pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
      year = add_level(N = (years*2), nest = FALSE),
      obs = cross_levels(
        by = join(pixels, year),
        post = ifelse(year > years, 1, 0),
        v_it = rnorm(N, 0, std_v),
        ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it
      )
    )
    
    #generate random 
    error_table <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(error_table, by = "property") %>%
      mutate(ystar = ystar + p_err) %>%
      mutate(y = (ystar > 0)*1 )
    
    #need to determine which year deforestation occurred
    year_df <- panels %>%
      dplyr::select(pixels, year, y) %>%
      dcast(pixels ~ year , value.var = "y")
    
    rownames(year_df) <- year_df$pixels
    
    year_df <- year_df %>%
      dplyr::select(- pixels)
    
    #creating variable for the year a pixel is deforested
    not_defor <- rowSums(year_df)<1 *1
    defor_year <- max.col(year_df, ties.method = "first") 
    defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
    defor_df <- tibble::rownames_to_column(defor_df)
    names(defor_df)[1] <- paste("pixels")
    
    panels <- defor_df %>%
      dplyr::select(pixels, defor_year) %>%
      inner_join(panels, by = "pixels")
    
    
    
    cols.num <- c("pixels", "grid", "property", "county", "year")
    panels[cols.num] <- sapply(panels[cols.num],as.numeric)
    
    panels <- panels %>%
      mutate(indic = year - defor_year) %>%
      mutate(defor = ifelse(indic > 0, 1, y))%>%
      mutate(y_it = ifelse(indic > 0, NA, y))
    
    
    DID <- lm(y_it ~  post*treat, 
               data   = panels
    )
    
    #calculating bias from each aggregation method
    coeffmatrix[i,1] <- DID$coefficients[4] - ATT
    
    
    
    #calculating standard errors and whether att is within CI
    cluster_county_pix    <- sqrt(diag(vcovCR(DID, panels$county, type="CR0")))[4]
    covermat[i,1] <- between(ATT, DID$coefficients[4] - 1.96 * cluster_county_pix, DID$coefficients[4] + 1.96 * cluster_county_pix)*1

    print(i)
    toc()
  }
  
  
  outputs = list("covermat"= covermat)
  return(outputs)
  
  #end function  
}  
