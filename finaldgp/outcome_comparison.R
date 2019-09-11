library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(tictoc)
source('grid_scapegen.R')

#begin function
outcome_comparison <- function(n, nobs, years, b0, b1, b2, b3, std_a = 0.1, std_v = 0.25, cellsize){
  
  gridscape = grid_scapegen(nobs, cellsize)
  pixloc_df = gridscape$pixloc_df
  gridcoords = gridscape$gridcoords
  #N_treat = gridscape$N_treat    
  ATT <- pnorm(b0+b1+b2+b3, 0, (std_a^2+std_v^2)^.5) - pnorm(b0+b1+b2, 0, (std_a^2+std_v^2)^.5)
  
  pixloc <- pixloc_df[order(pixloc_df$pixels),]
  
  coeffmatrix <- matrix(nrow = n, ncol = 2)
  
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
    
    
    
    
    
    gridlevel_df <- gridlevel_df[order(gridlevel_df$grid, gridlevel_df$year),]
    gridlevel_df <- slide(gridlevel_df, Var = "defor", GroupVar = "grid", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    gridlevel_df$forshare <- 1 -  gridlevel_df$defor
    gridlevel_df$forsharelag <- 1 -  gridlevel_df$deforlag
    
    #### create baseline defor and forshare vars 
    #county level
    year1 <- subset(gridlevel_df, year== 1)
    colnames(year1)[colnames(year1)=="defor"] <- "defor0"
    year1 <- subset(year1, select = c(grid, defor0))
    #st_geometry(year1) <- NULL
    gridlevel_df <-  merge(gridlevel_df, year1, by = "grid")
    
    gridlevel_df$forshare0 <- 1 -  gridlevel_df$defor0
    
    #generate outcome var
    gridlevel_df$deforrate1 <- ((gridlevel_df$forsharelag- gridlevel_df$forshare) / gridlevel_df$forsharelag)
    
    gridlevel_df$deforrate2 <- gridlevel_df$forshare / gridlevel_df$forshare0
    
    #remove any infinite values
    gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
    gridlevel_df <- 
      gridlevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # run two-way fixed effects with outcome 1 
    coeffmatrix[i,1] <- plm(deforrate1 ~  post*treat, 
                            data   = gridlevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #grid and year fixed effects
                            index  = c("grid", "year")
    )$coefficients - ATT
    
    coeffmatrix[i,2] <- plm(deforrate2 ~  post*treat, 
                            data   = gridlevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #grid and year fixed effects
                            index  = c("grid", "year")
    )$coefficients - ATT
    
    print(i)
    toc()
  }
  
  
  coeff_bias <- as.data.frame(coeffmatrix)
  actual <- rep(ATT, times = n)
  names(coeff_bias)[1] <- paste("outcome1")
  names(coeff_bias)[2] <- paste("outcome2")
  suppressWarnings(cbias <- melt(coeff_bias, value.name = "bias"))
  
  
  plot <- ggplot(data = cbias, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("outcome1", "outcome2"), labels=c("outcome 2", "outcome 2"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Outcome 1 Mean:", round(mean(coeff_bias$outcome1), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$outcome1), digits = 4), "\n", 
                                    "Outcome 2 Mean:", round(mean(coeff_bias$outcome2), digits = 4),
                                    ", RMSE:", round(rmse(actual, coeff_bias$outcome2), digits = 4) ) 
    )
  
  
  
  outputs = list("plot" = plot, "biases" = coeff_bias)
  return(outputs)
  
  #end function  
}  
