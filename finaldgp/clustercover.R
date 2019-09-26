#this function uses the property_scape gen fcn to generate a landscape with properties. Then we introduce property level perturbations and compare estimates when the data is aggregateed to the grid vs. property level


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
clustercover <- function(n, nobs, years, b0, b1, b2, b3, std_a = 0.1, std_v = 0.25, std_p = .1, cellsize, ppoints, cpoints){
  
  countyscape = county_scapegen(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  ATT <- pnorm(b0+b1+b2+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  DID_estimand <- (pnorm(b0+b1+b2+b3, 0, (std_a^2+std_v^2+std_p^2)^.5) - pnorm(b0+b1, 0, (std_a^2+std_v^2+std_p^2)^.5)
                   - (pnorm(b0+b2, 0, (std_a^2+std_v^2+std_p^2)^.5) - pnorm(b0, 0, (std_a^2+std_v^2+std_p^2)^.5)) )
  
  pixloc <- pixloc_df[order(pixloc_df$pixels),]
  covermat <- matrix(nrow = n, ncol = 6)
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  
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
        ystar = b0 + b1*treat + b2*post + b3*treat*post + a_i + v_it
      )
    )
    
    #generate random 
    error_table <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    panels <- inner_join(pixloc, panels, by = c("pixels", "treat"))
    panels <- inner_join(panels, error_table, by = c("property"))
    
    panels$ystar <- panels$ystar + panels$p_err
    panels$y <- (panels$ystar > 0)*1
    
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
    panels <- inner_join(defor_df, panels, by = "pixels")
    
    
    # creating three outcome variables for each possible situation
    ### y: allows the outcome to switch between 0 and 1 across years
    ### defor: outcome is set to 1 in each year after the pixel is deforested
    #sapply(pixloc_df, class)
    cols.num <- c("pixels", "grid", "property", "county", "year")
    panels[cols.num] <- sapply(panels[cols.num],as.numeric)
    
    panels$indic <- (panels$year - panels$defor_year)
    panels$defor <- ifelse(panels$indic > 0 , 1, panels$y)
    #panels <- subset(panels, select = -c(indic))
    
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
    
    #calculating bias from each aggregation method
    coeffmatrix[i,1] <- DID1$coefficients - ATT
    coeffmatrix[i,2] <- DID2$coefficients - ATT
    coeffmatrix[i,3] <- DID3$coefficients - ATT
    
    cluster_se1    <- sqrt(diag(vcovHC(DID1, type = "HC0", cluster = "group")))
    cluster_se2    <- sqrt(diag(vcovHC(DID2, type = "HC0", cluster = "group")))
    cluster_se3    <- sqrt(diag(vcovHC(DID3, type = "HC0", cluster = "group")))
    covermat[i,1] <- between(ATT, DID1$coefficients - 1.96 * cluster_se1, DID1$coefficients + 1.96 * cluster_se1)*1
    covermat[i,2] <- between(ATT, DID2$coefficients - 1.96 * cluster_se2, DID2$coefficients + 1.96 * cluster_se2)*1
    covermat[i,3] <- between(ATT, DID3$coefficients - 1.96 * cluster_se3, DID3$coefficients + 1.96 * cluster_se3)*1
    
    se1 <- sqrt(diag(vcovHC(DID1, type = "HC0")))
    se2 <- sqrt(diag(vcovHC(DID2, type = "HC0")))
    se3 <- sqrt(diag(vcovHC(DID3, type = "HC0")))
    covermat[i,4] <- between(ATT, DID1$coefficients - 1.96 * se1, DID1$coefficients + 1.96 * se1)*1
    covermat[i,5] <- between(ATT, DID2$coefficients - 1.96 * se2, DID2$coefficients + 1.96 * se2)*1
    covermat[i,6] <- between(ATT, DID3$coefficients - 1.96 * se3, DID3$coefficients + 1.96 * se3)*1
    
    
    print(i)
    toc()
  }
  
  coeff_bias <- as.data.frame(coeffmatrix)
  actual <- rep(ATT, times = n)
  names(coeff_bias)[1] <- paste("grid")
  names(coeff_bias)[2] <- paste("property")
  names(coeff_bias)[3] <- paste("county")
  suppressWarnings(cbias <- melt(coeff_bias, value.name = "bias"))
  
  grid_clustercover <- mean(covermat[,1]) 
  prop_clustercover <- mean(covermat[,2])
  county_clustercover <- mean(covermat[,3])  
  
  grid_cover <- mean(covermat[,4]) 
  prop_cover <- mean(covermat[,5])
  county_cover <- mean(covermat[,6])  
  
  
  plot <- ggplot(data = cbias, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("grid", "property", "county"), labels=c("aggregated to grids", "aggregated to properties", "aggregated to counties"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #geom_vline(aes(xintercept= (DID_estimand - ATT), color="DID estimand - ATT"), linetype="dashed")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Mean grid:", round(mean(coeff_bias$grid), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$grid), digits = 4), "\n", 
                                    "Mean property:", round(mean(coeff_bias$property), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$property), digits = 4),"\n",
                                    "Mean county:", round(mean(coeff_bias$county), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$county), digits = 4) ) 
    )
  
  outputs = list("plot" = plot, "biases" = coeff_bias, "grid_cover" = grid_cover, "prop_cover" = prop_cover, "county_cover" = county_cover, "grid_clustercover" = grid_clustercover, "prop_clustercover" = prop_clustercover, "county_clustercover" = county_clustercover)
  return(outputs)
  
  #end function  
}  
