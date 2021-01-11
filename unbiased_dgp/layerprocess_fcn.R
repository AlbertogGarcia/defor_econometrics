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
library(fixest)
source('full_landscape.R')

#begin function
layerprocess_fcn <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, std_p = 0.0, std_c = 0.0, cellsize, ppoints, cpoints){
  
  countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2+std_c^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2+std_c^2)^.5)
  
  pixloc <- pixloc_df
  
  covermat <- matrix(nrow = n, ncol = 6)
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  
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
    errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    errortable_county <- data.frame(county = as.character(unique(pixloc$county)), c_err = rnorm(length(unique(pixloc$county)), 0, std_c))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(errortable_property, by = "property") %>%
      inner_join(errortable_county, by = "county") %>%
      mutate(ystar = ystar + p_err + c_err) %>%
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
    
    # aggregate up to county in each year 
    suppressWarnings(
      gridlevel_df <-  aggregate(panels, by = list(panels$year, panels$grid), FUN = mean, drop = TRUE)[c("grid", "treat", "post", "year","defor", "garea")]
    )
    
    gridlevel_df <- gridlevel_df[order(gridlevel_df$grid, gridlevel_df$year),]
    gridlevel_df <- slide(gridlevel_df, Var = "defor", GroupVar = "grid", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    gridlevel_df$forshare <- 1 -  gridlevel_df$defor
    gridlevel_df$forsharelag <- 1 -  gridlevel_df$deforlag
    
    #generate outcome var
    gridlevel_df$deforrate <- ((gridlevel_df$forsharelag- gridlevel_df$forshare) / gridlevel_df$forsharelag)
    #remove any infinite values
    #gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
    gridlevel_df <- 
      gridlevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregate up to property in each year 
    suppressWarnings(
      proplevel_df <-  aggregate(panels, by = list(panels$property, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("property", "treat", "post", "year","defor", "parea")]
    )
    
    
    proplevel_df <- proplevel_df[order(proplevel_df$property, proplevel_df$year),]
    proplevel_df <- slide(proplevel_df, Var = "defor", GroupVar = "property", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    proplevel_df$forshare <- 1 -  proplevel_df$defor
    proplevel_df$forsharelag <- 1 -  proplevel_df$deforlag
    
    #generate outcome var
    proplevel_df$deforrate <- ((proplevel_df$forsharelag- proplevel_df$forshare) / proplevel_df$forsharelag)
    #remove any infinite values
    
    proplevel_df <- 
      proplevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # aggregate up to county in each year 
    suppressWarnings(
      countylevel_df <-  aggregate(panels, by = list(panels$county, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("county", "property", "treat", "post", "year","defor", "carea")]
    )
    
    
    countylevel_df <- countylevel_df[order(countylevel_df$county, countylevel_df$year),]
    countylevel_df <- slide(countylevel_df, Var = "defor", GroupVar = "county", NewVar = "deforlag",
                            slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    countylevel_df$forshare <- 1 -  countylevel_df$defor
    countylevel_df$forsharelag <- 1 -  countylevel_df$deforlag
    
    #generate outcome var
    countylevel_df$deforrate <- ((countylevel_df$forsharelag- countylevel_df$forshare) / countylevel_df$forsharelag)
    #remove any infinite values
    
    countylevel_df <- 
      countylevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    
    # DID1 <- plm(deforrate ~  post*treat,
    #             data   = gridlevel_df,
    #             #weights=(gridlevel_df$garea),
    #             method = "within", #fixed effects model
    #             effect = "twoway", #property and year fixed effects
    #             index  = c("grid", "year")
    # )
    
    DID1 <- feols(deforrate ~  post*treat|year+grid, data = gridlevel_df)
    
    # run two-way fixed effects with outcome 1 
    # DID2 <- plm(deforrate ~  post*treat, 
    #             data   = proplevel_df, 
    #             #weights=(proplevel_df$parea),
    #             method = "within", #fixed effects model
    #             effect = "twoway", #property and year fixed effects
    #             index  = c("property", "year")
    # )
    
    DID2 <- feols(deforrate ~  post*treat|year+property, data = proplevel_df)
    
    
    # DID3 <- plm(deforrate ~  post*treat, 
    #             data   = countylevel_df, 
    #             #weights=(countylevel_df$carea),
    #             method = "within", #fixed effects model
    #             effect = "twoway", #county and year fixed effects
    #             index  = c("county", "year")
    # )
    
    DID3 <- feols(deforrate ~  post*treat|year+county, data = countylevel_df)
    
    
    # DID4 <- lm(y_it ~  post*treat, 
    #            data   = panels
    # )
    
    DID4 <- feols(y_it ~  post*treat, data = panels)
    
    
    #calculating bias from each aggregation method
    coeffmatrix[i,1] <- DID1$coefficients - ATT
    coeffmatrix[i,2] <- DID2$coefficients - ATT
    coeffmatrix[i,3] <- DID3$coefficients - ATT
    coeffmatrix[i,4] <- DID4$coefficients[4] - ATT
    
    #clustering at group level for aggregated analyses
    cluster_se1    <- summary(DID1, cluster = ~grid)$se[1]
    cluster_se2    <- summary(DID2, cluster = ~property)$se[1]
    cluster_se3    <- summary(DID3, cluster = ~county)$se[1]
    
    # clustered pixel level standard errors
    cluster_prop_pix  <- summary(DID4, cluster = ~property)$se[4]
    cluster_grid_pix <- summary(DID4, cluster = ~grid)$se[4]
    cluster_county_pix <- summary(DID4, cluster = ~county)$se[4]
    
    #whether att is within CI
    
    # coverage with clustered standard errors at group level
    covermat[i,1] <- between(ATT, DID1$coefficients - 1.96 * cluster_se1, DID1$coefficients + 1.96 * cluster_se1)*1
    covermat[i,2] <- between(ATT, DID2$coefficients - 1.96 * cluster_se2, DID2$coefficients + 1.96 * cluster_se2)*1
    covermat[i,3] <- between(ATT, DID3$coefficients - 1.96 * cluster_se3, DID3$coefficients + 1.96 * cluster_se3)*1
    
    # coverage with pixel level analyses clustered at different levels
    covermat[i,4] <- between(ATT, DID4$coefficients[4] - 1.96 * cluster_prop_pix, DID4$coefficients[4] + 1.96 * cluster_prop_pix)*1
    covermat[i,5] <- between(ATT, DID4$coefficients[4] - 1.96 * cluster_grid_pix, DID4$coefficients[4] + 1.96 * cluster_grid_pix)*1
    covermat[i,6] <- between(ATT, DID4$coefficients[4] - 1.96 * cluster_county_pix, DID4$coefficients[4] + 1.96 * cluster_county_pix)*1
    
    print(i)
    toc()
  }
  
  coeff_bias <- as.data.frame(coeffmatrix)
  actual <- rep(0, times = n)
  names(coeff_bias)[1] <- paste("grid")
  names(coeff_bias)[2] <- paste("property")
  names(coeff_bias)[3] <- paste("county")
  names(coeff_bias)[4] <- paste("pixel")
  
  
  aggregation <- c('grid', 
                   'property',
                   'county',
                   'pixel',
                   'pixel',
                   'pixel'
  )
  
  std_error <- c('clustered at grid', 
                 'clustered at property',
                 'clustered at county',
                 'clustered at property',
                 'clustered at grid',
                 'clustered at county'
  )
  
  coverages_df <- data.frame(aggregation, std_error, coverage = colMeans(covermat)) 
  
  
  outputs = list("plot" = plot, "biases" = coeff_bias, "coverages_df" = coverages_df)
  return(outputs)
  
  #end function  
}  
