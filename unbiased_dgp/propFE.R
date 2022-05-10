#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level

library(ggplot2)
library(clubSandwich)
library(reshape2)
library(matrixStats)
library(ggplot2)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tidyverse)
library(tictoc)
library(fixest)
library(here)
library(DeclareDesign)
library(msm)
library(survival)

source(here::here('unbiased_dgp', 'multigrid_landscape.R'))
source(here::here('unbiased_dgp', 'proptreat_landscape.R'))

#begin function
propFE <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.5, std_p = 0.0, cellsize_list, ppoints, cpoints, proptreatassign = FALSE){
  
  if (proptreatassign) {
    countyscape = proptreat_landscape(nobs, cellsize_list, ppoints, cpoints)
    pixloc_df = countyscape$pixloc_df
  } else {
    countyscape = multigrid_landscape(nobs, cellsize_list, ppoints, cpoints)
    pixloc_df = countyscape$pixloc_df
  }
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  pixloc <- pixloc_df
  
  n_mod = 3
  
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'b3'= rep(b3, summ_row), 
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row), "years" = years,
                             'iteration' = rep(NA, summ_row), 
                             'pixel'=rep(NA, summ_row), 'property fe'=rep(NA, summ_row), 'treatment var' = rep(NA, summ_row), 'se_property'=rep(NA, summ_row),
                             'bias'=rep(NA, summ_row), 'cover'=rep(NA, summ_row),'notes'=rep(NA, summ_row),
                             stringsAsFactors=FALSE)
  
  options(dplyr.summarise.inform = FALSE)
  print(min(pixloc$parea))
  
  for(i in 1:n){
    tic("loop")
    Nobs <- length(pixloc$treat)  
    panels <- fabricate(
      pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
      year = add_level(N = (years*2), nest = FALSE),
      obs = cross_levels(
        by = join(pixels, year),
        post = ifelse(as.numeric(year) > years, 1, 0),
        v_it = rnorm(N, 0, std_v),
        ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
        ystar_cf = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
      )
    )
    
    #generate random 
    errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(errortable_property, by = "property") %>%
      #inner_join(errortable_county, by = "county") %>%
      mutate(year = as.numeric(year),
             ystar = ystar + p_err,
             ystar_cf = ystar_cf + p_err,
             y = (ystar > 0)*1 ,
             y_cf = (ystar_cf > 0)*1 )
    
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
      inner_join(panels, by = "pixels") %>%
      mutate_at(vars(pixels, county, year, property, defor_year), as.numeric)
    
    panels <- panels %>%
      mutate(indic = year - defor_year,
             defor = ifelse(indic > 0, 1, y),
             y_it = ifelse(indic > 0, NA, y))
    
    panels <- panels %>%
      group_by(property, year)%>%
      mutate(ptreat = mean(treat))
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ##### RUNNING MODELS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ### TWFE regressions with aggregated fixed effects
    
    ## binary treatment
    propfe1 <- suppressMessages(feols(y_it ~  post:treat|year + property, data = panels))
    
    propfe1_se    <- tail(summary(propfe1, cluster = ~property)$se, n=1)
    
    propfe1_cover <- dplyr::between(ATT, tail(propfe1$coefficients, n=1) - 1.96 * propfe1_se, tail(propfe1$coefficients, n=1) + 1.96 * propfe1_se)*1
    
    
    ## binary treatment + treat
    propfe2 <- suppressMessages(feols(y_it ~  post*treat|year + property, data = panels))
    
    propfe2_se    <- tail(summary(propfe2, cluster = ~property)$se, n=1)
    
    propfe2_cover <- dplyr::between(ATT, tail(propfe2$coefficients, n=1) - 1.96 * propfe2_se, tail(propfe2$coefficients, n=1) + 1.96 * propfe2_se)*1
    
   
    
    ## continuous treatment
    propfe3 <- suppressMessages(feols(y_it ~  post:ptreat|year + property, data = panels))
    
    propfe3_se    <- tail(summary(propfe3, cluster = ~property)$se, n=1)
    
    propfe3_cover <- dplyr::between(ATT, tail(propfe3$coefficients, n=1) - 1.96 * propfe3_se, tail(propfe3$coefficients, n=1) + 1.96 * propfe3_se)*1
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ### filling in summary long dataframe output
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="notes")
    
    # Pixel fixed effects
    summary_long[i,c(firstcol:lastcol)] <- c(
      i,
      tail(propfe1$coefficients, n=1) - ATT,
      propfe1_cover,
      "binary"
    )
    
    # traditional DID
    summary_long[i+n,c(firstcol:lastcol)] <- c(
      i,
      tail(propfe2$coefficients, n=1) - ATT,
      propfe2_cover,
      "binary + treat"
    )
    
    # property uoa
    summary_long[i+n*2,c(firstcol:lastcol)] <- c(
      i,
      tail(propfe3$coefficients, n=1) - ATT,
      propfe3_cover,
      "continuous"
    )
    
    
    print(i)
    toc()
  }
  
  outputs = list(
    "summary_long" = summary_long)
  
  return(outputs)
  
  #end function  
}  