#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level

library(ggplot2)
library(clubSandwich)
library(reshape2)
library(matrixStats)
library(ggplot2)
# library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tidyverse)
library(tictoc)
library(fixest)
library(here)
library(DeclareDesign)
source(here::here('unbiased_dgp', 'multigrid_landscape.R'))

#begin function
multiple_gridsize_clean <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, std_p = 0.0, cellsize_list, ppoints, cpoints, rm.selection = FALSE){
  
  countyscape = multigrid_landscape(nobs, cellsize_list, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  
  pixloc <- pixloc_df
  
  prop_area <- pixloc %>%
    group_by(property) %>%
    summarise(parea = mean(parea))
  
  grid_covermat <- matrix(nrow = n, ncol = length(cellsize_list))
  grid_coeffmatrix <- grid_covermat
  pix_coeffmatrix <- grid_covermat
  pix_covermat <- grid_covermat
  
  n_mod = length(cellsize_list)*2
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row),  
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row),
                             'grid fe'=rep(1, summ_row), 'parea' = rep(mean(prop_area$parea), summ_row),
                             'pixel'=rep(NA, summ_row), 'grid'=rep(NA, summ_row),
                             'iteration' = rep(NA, summ_row), 'avg_garea' = rep(NA, summ_row),
                             'bias'=rep(NA, summ_row), 'cover'=rep(NA, summ_row),
                             stringsAsFactors=FALSE)
  
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
        ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
        ystar_cf = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
      )
    )
    
    #generate random 
    errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    #errortable_county <- data.frame(county = as.character(unique(pixloc$county)), c_err = rnorm(length(unique(pixloc$county)), 0, std_c))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(errortable_property, by = "property") %>%
      #inner_join(errortable_county, by = "county") %>%
      mutate(ystar = ystar + p_err,
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
      mutate(indic = year - defor_year) %>%
      mutate(defor = ifelse(indic > 0, 1, y))%>%
      mutate(y_it = ifelse(indic > 0, NA, y))
    
    treat_post <- subset(panels, treat==1)
    
    control_post <- subset(panels, treat==0 )
    
    y_treat = mean(treat_post$y)
    y0_treat = mean(treat_post$y_it, na.rm = TRUE)
    y_control = mean(control_post$y)
    y0_control = mean(control_post$y_it, na.rm = TRUE)
    
    
    sel_bias = y0_treat - y_treat - (y0_control - y_control)
    
    ATT = mean( subset(panels, treat==1&post==1)$y)-mean( subset(panels, treat==1&post==1)$y_cf)
    
    if(rm.selection){
      ATT = ATT+sel_bias  
    }
    
    for(k in 1:length(cellsize_list)){
      
      this_panel <- panels %>%
        rename(this_grid = paste0("grid_", cellsize_list[[k]]),
               this_garea = paste0("garea_", cellsize_list[[k]]))%>%
        select(this_grid, treat, post, year, defor, this_garea, y_it)
      
      # aggregate up to county in each year 
      this_grid_df <- as.data.frame(this_panel) %>%
        dplyr::group_by(this_grid, year, post) %>%
        dplyr::summarise(this_garea = mean(this_garea),
                         defor = mean(defor),
                         treat = mean(treat))%>%
        ungroup()
      
      # suppressWarnings(
      #   this_grid_df <-  aggregate(this_panel, by = list(this_panel$year, this_panel$this_grid), FUN = mean, drop = TRUE)
      #   )
      
      this_grid_df <- this_grid_df[order(this_grid_df$this_grid, this_grid_df$year),]
      this_grid_df <- slide(this_grid_df, Var = "defor", GroupVar = "this_grid", NewVar = "deforlag",
                            slideBy = -1, reminder = FALSE)
      
      ##### creating forested share variable #####
      this_grid_df$forshare <- 1 -  this_grid_df$defor
      this_grid_df$forsharelag <- 1 -  this_grid_df$deforlag
      
      #generate outcome var
      this_grid_df$deforrate <- ((this_grid_df$forsharelag- this_grid_df$forshare) / this_grid_df$forsharelag)
      #remove any infinite values
      #gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
      this_grid_df <- 
        this_grid_df %>% 
        filter_all(all_vars(!is.infinite(.)))
      
      
      # aggregated units of analysis
      DID <- feols(deforrate ~  post*treat|year+this_grid, data = this_grid_df)
      
      # bias
      grid_coeffmatrix[i,k] <- DID$coefficients - ATT
      
      #clustering at grid level
      clse    <- tail(summary(DID, cluster = ~this_grid)$se, n=1)
      
      grid_covermat[i,k] <- dplyr::between(ATT, tail(DID$coefficients, n=1) - 1.96 * clse, tail(DID$coefficients, n=1) + 1.96 * clse)*1
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      # aggregated units of analysis
      pix_twfe <- feols(y_it ~  post*treat|year+this_grid, data = this_panel)
      
      # bias
      pix_coeffmatrix[i,k] <- tail(pix_twfe$coefficients, n=1) - ATT
      
      #clustering at grid level
      clse    <- tail(summary(pix_twfe, cluster = ~this_grid)$se, n=1)
      
      pix_covermat[i,k] <- dplyr::between(ATT, tail(pix_twfe$coefficients, n=1) - 1.96 * clse, tail(pix_twfe$coefficients, n=1) + 1.96 * clse)*1
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      firstcol = which(colnames(summary_long)=="pixel")
      lastcol = which(colnames(summary_long)=="cover")
      
      summary_long[i+((k-1)*n),c(firstcol:lastcol)] <- c(
        0,1,
        i,
        mean(this_grid_df$this_garea), 
        grid_coeffmatrix[i,k],
        grid_covermat[i,k]
      )
      
      summary_long[i+((k-1)*n)+(summ_row/2),c(firstcol:lastcol)] <- c(
        1,0,
        i,
        mean(this_grid_df$this_garea), 
        pix_coeffmatrix[i,k],
        pix_covermat[i,k]
      )
      
    }
    
    
    print(i)
    toc()
  }
  
  
  summary_long <- summary_long %>%
    select(iteration, everything())
  
  
  outputs = list(
    "summary_long" = summary_long,
    "pixloc_df" = pixloc##,
   # "selection_bias" = selection_bias
   )
  return(outputs)
  
  #end function  
}  
