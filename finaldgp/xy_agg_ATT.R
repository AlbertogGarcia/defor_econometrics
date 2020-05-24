library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tictoc)
source('county_scapegen.R')

#begin function
xy_agg_ATT <- function(n, nobs, years, min_ATT, max_ATT, base_0, base_1, trend, std_a = 0.1, std_v = 0.25, std_p = .1, cellsize, ppoints, cpoints){
  
  ATT_vals <- seq(from = min_ATT, to = max_ATT, by = .01)
  list <- as.list(ATT_vals)
  
  #preallocate n x ___ matrix
  gridmatrix <- matrix(nrow = n, ncol = length(list))
  countymatrix <- matrix(nrow = n, ncol = length(list))
  propertymatrix <- matrix(nrow = n, ncol = length(list))
  
  std_avp = (std_a^2 + std_v^2 + std_p)^.5
  b0 = qnorm(base_0, mean = 0, sd = std_avp)
  b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
  b2 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
  
  countyscape = county_scapegen(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  
  pixloc <- pixloc_df[order(pixloc_df$pixels),]
  
  for(k in list){
    tic("loop")
    b3 = qnorm( pnorm(b0+b1+b2, mean = 0, sd = std_avp) + k , mean = 0, sd = std_avp) - (b0 + b1 + b2)
    ATT = pnorm(b0+b1+b2+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2, 0, (std_a^2+std_v^2 + std_p^2)^.5)
    
    for(i in 1:n){
      Nobs <- length(pixloc$treat)  
      panels <- fabricate(
        pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
        year = add_level(N = (years*2), nest = FALSE),
        obs = cross_levels(
          by = join(pixels, year),
          post = ifelse(year > years, 1, 0),
          v_it = rnorm(N, 0, std_v),
          ystar = b0 + b1*treat + b2*post + b3*treat*post + a_i + v_it
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
        select(pixels, year, y) %>%
        dcast(pixels ~ year , value.var = "y")
      
      rownames(year_df) <- year_df$pixels
      
      year_df <- year_df %>%
        select(- pixels)
      
      #creating variable for the year a pixel is deforested
      not_defor <- rowSums(year_df)<1 *1
      defor_year <- max.col(year_df, ties.method = "first") 
      defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
      defor_df <- tibble::rownames_to_column(defor_df)
      names(defor_df)[1] <- paste("pixels")
      
      panels <- defor_df %>%
        select(pixels, defor_year) %>%
        inner_join(panels, by = "pixels")
      
      
      
      cols.num <- c("pixels", "grid", "property", "county", "year")
      panels[cols.num] <- sapply(panels[cols.num],as.numeric)
      
      panels <- panels %>%
        mutate(indic = year - defor_year) %>%
        mutate(defor = ifelse(indic > 0, 1, y)
        )
      
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
        countylevel_df <-  aggregate(panels, by = list(panels$county, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("county", "treat", "post", "year","defor", "carea")]
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
      
      
      
      
      
      
      DID1 <- plm(deforrate ~  post*treat, 
                  data   = gridlevel_df, 
                  #weights=(gridlevel_df$garea),
                  method = "within", #fixed effects model
                  effect = "twoway", #property and year fixed effects
                  index  = c("grid", "year")
      )
      
      # run two-way fixed effects with outcome 1 
      DID2 <- plm(deforrate ~  post*treat, 
                  data   = proplevel_df, 
                  #weights=(proplevel_df$parea),
                  method = "within", #fixed effects model
                  effect = "twoway", #property and year fixed effects
                  index  = c("property", "year")
      )
      
      DID3 <- plm(deforrate ~  post*treat, 
                  data   = countylevel_df, 
                  #weights=(countylevel_df$carea),
                  method = "within", #fixed effects model
                  effect = "twoway", #county and year fixed effects
                  index  = c("county", "year")
      )
      
      place <- which(list == k)
      
      #calculating bias from each aggregation method
      gridmatrix[i,place] <- DID1$coefficients - ATT
      countymatrix[i,place] <- DID3$coefficients - ATT
      propertymatrix[i,place] <- DID2$coefficients - ATT
      
    }
    print(k)
    toc()
  }
  
  gridbias_df <- as.data.frame(cbind(colMeans(gridmatrix), ATT_vals)) 
  
  countybias_df <- as.data.frame(cbind(colMeans(countymatrix), ATT_vals)) 
  
  propbias_df <- as.data.frame(cbind(colMeans(propertymatrix), ATT_vals)) 
  
  
  plot <- ggplot() + 
    geom_line(data = gridbias_df, aes(x = ATT_vals, y = V1, colour="aggregated to grid"), color = "purple", size =1) +
    geom_line(data = countybias_df, aes(x = ATT_vals, y = V1, colour="aggregated to county"), color = "green", size =1) +
    geom_line(data = propbias_df, aes(x = ATT_vals, y = V1, colour="aggregated to property"), color = "blue", size =1) +
    labs(x = "ATT", y = "Bias", caption = paste("Bias as a function of ATT depends on aggregation method")) + 
    geom_hline(yintercept = 0, linetype = "dashed")
  
  
  outputs = list("gridbias_df" = gridbias_df, "countybias_df" = countybias_df, "propbias_df" = propbias_df, "plot" = plot)
  
  return(outputs)
  
}
    
      
    
