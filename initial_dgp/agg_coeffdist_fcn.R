# function designed to give mean and variance of coefficients when aggregated to county or property level

library(clubSandwich)
library(matrixStats)
library(DataCombine)
library(plm)
library(IDPmisc)
library(tictoc)
#begin function
agg_coeffdist_fcn <- function(n, nobs, years, psize, cellsize){
  
  #preallocate n x _ matrix
  coeffmatrix <- matrix(nrow = n, ncol = 14)
  
  for(i in 1:n){
    
    tic("loop")
    # call defor_sim function to simulate dataframe, returned as defor_df  
    prop_deforsim(nobs, years, psize, cellsize)  
    #prop_deforsim(10000,3,1200,25)  
    ############### gather area for each property and county ##################
    
    
    

    # lagging deforestation and forest share variable

    countypert_df <- countypert_df[order(countypert_df$countyid, countypert_df$year),]
    countypert_df <- slide(countypert_df, Var = "defor", GroupVar = "countyid", NewVar = "deforlag",
                           slideBy = -1, reminder = FALSE)
    #property level
    proppert_df <- proppert_df[order(proppert_df$property, proppert_df$year),]
    proppert_df <- slide(proppert_df, Var = "defor", GroupVar = "property", NewVar = "deforlag",
                         slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    countypert_df$forshare <- 1 -  countypert_df$defor
    proppert_df$forshare <- 1 - proppert_df$defor
    countypert_df$forsharelag <- 1 -  countypert_df$deforlag
    proppert_df$forsharelag <- 1 - proppert_df$deforlag
    
    #### create baseline defor and forshare vars 
    #county level
    year1 <- subset(countypert_df, year== 1)
    colnames(year1)[colnames(year1)=="defor"] <- "defor0"
    year1 <- subset(year1, select = c(countyid, defor0))
    #st_geometry(year1) <- NULL
    countypert_df <-  merge(countypert_df, year1, by = "countyid")
    #property level
    year1 <- subset(proppert_df, year== 1)
    colnames(year1)[colnames(year1)=="defor"] <- "defor0"
    year1 <- subset(year1, select = c(property, defor0))
    #st_geometry(year1) <- NULL
    proppert_df <-  merge(proppert_df, year1, by = "property")
    countypert_df$forshare0 <- 1 -  countypert_df$defor0
    proppert_df$forshare0 <- 1 - proppert_df$defor0
    
    
    ##### constructing various outcome variables #####
    
    
    countypert_df$deforrate1 <- countypert_df$defor / countypert_df$defor0
    proppert_df$deforrate1 <- proppert_df$defor / proppert_df$defor0
    
    countypert_df$deforrate2 <- countypert_df$defor / countypert_df$deforlag
    proppert_df$deforrate2 <- proppert_df$defor / proppert_df$deforlag
    
    countypert_df$deforrate3 <- countypert_df$defor / countypert_df$forshare0
    proppert_df$deforrate3 <- proppert_df$defor / proppert_df$forshare0
    
    countypert_df$deforrate4 <- countypert_df$defor / countypert_df$forsharelag
    proppert_df$deforrate4 <- proppert_df$defor / proppert_df$forsharelag
    
    countypert_df$deforrate6 <- countypert_df$forshare / countypert_df$forshare0
    proppert_df$deforrate6 <- proppert_df$forshare / proppert_df$forshare0
    
    countypert_df$deforrate7 <- countypert_df$forshare / countypert_df$forsharelag
    proppert_df$deforrate7 <- proppert_df$forshare / proppert_df$forsharelag
    
    #remove infinite values

    proppert_df <- subset(proppert_df, select = -c(geometry))
    proppert_df <- 
    proppert_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    countypert_df <- subset(countypert_df, select = -c(geometry))
    countypert_df <- 
      countypert_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    # run two-way fixed effects with propertylevel data
    coeffmatrix[i,1] <- plm(deforrate1 ~  post*treat, 
                            data   = proppert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #property and year fixed effects
                            index  = c("property", "year")
    )$coefficients
    
    # run two-way fixed effects with county level data
    coeffmatrix[i,2] <- plm(deforrate1 ~  post*treat, 
                            data   = countypert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #county and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    
    # # run two-way fixed effects with propertylevel data
    # coeffmatrix[i,3] <- plm(deforrate2 ~  post*treat, 
    #                         data   = proppert_df, 
    #                         method = "within", #fixed effects model
    #                         effect = "twoway", #property and year fixed effects
    #                         index  = c("property", "year")
    # )$coefficients
    # 
    # # run two-way fixed effects with county level data
    # coeffmatrix[i,4] <- plm(deforrate2 ~  post*treat, 
    #                         data   = countypert_df, 
    #                         method = "within", #fixed effects model
    #                         effect = "twoway", #county and year fixed effects
    #                         index  = c("countyid", "year")
    # )$coefficients
    # 
    
    coeffmatrix[i,3]  <- lm(deforrate2 ~  post*treat, 
                            data = countypert_df
    )$coefficients[4]
    
    coeffmatrix[i,4]  <- lm(deforrate2 ~  post*treat, 
                            data = proppert_df
    )$coefficients[4]
    
    
    
    # run two-way fixed effects with propertylevel data
    coeffmatrix[i,5] <- plm(deforrate3 ~  post*treat, 
                            data   = proppert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #property and year fixed effects
                            index  = c("property", "year")
    )$coefficients
    
    # run two-way fixed effects with county level data
    coeffmatrix[i,6] <- plm(deforrate3 ~  post*treat, 
                            data   = countypert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #county and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients

 
      # run two-way fixed effects with propertylevel data
      # coeffmatrix[i,7] <- plm(deforrate4 ~  post*treat, 
      #                         data   = proppert_df, 
      #                         method = "within", #fixed effects model
      #                         effect = "twoway", #property and year fixed effects
      #                         index  = c("property", "year")
      # )$coefficients
    
    # run two-way fixed effects with county level data
    # coeffmatrix[i,8] <- plm(deforrate4 ~  post*treat, 
    #                         data   = countypert_df, 
    #                         method = "within", #fixed effects model
    #                         effect = "twoway", #county and year fixed effects
    #                         index  = c("countyid", "year")
    # )$coefficients
    
    coeffmatrix[i,7]  <- lm(deforrate4 ~  post*treat, 
                            data = countypert_df
    )$coefficients[4]
    
    coeffmatrix[i,8]  <- lm(deforrate4 ~  post*treat, 
                            data = proppert_df
    )$coefficients[4]
    
    ######################### Busch et al (2014): (f_it-1 -f_it)/f_it-1  ########################### 
    ### creating deforrate5 outcome variable for both dataframes
    
    proppert_df$deforrate5 <- ((proppert_df$forsharelag - proppert_df$forshare) / proppert_df$forsharelag)
    
    
    countypert_df$deforrate5 <- ((countypert_df$forsharelag- countypert_df$forshare) / countypert_df$forsharelag)
    
      
      # run two-way fixed effects with propertylevel data
      coeffmatrix[i,9] <- plm(deforrate5 ~  post*treat, 
                              data   = proppert_df, 
                              method = "within", #fixed effects model
                              effect = "twoway", #property and year fixed effects
                              index  = c("property", "year")
      )$coefficients
    
    # run two-way fixed effects with county level data
    coeffmatrix[i,10] <- plm(deforrate5 ~  post*treat, 
                            data   = countypert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #county and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    # run two-way fixed effects with propertylevel data
    coeffmatrix[i,11] <- plm(deforrate6 ~  post*treat, 
                            data   = proppert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #property and year fixed effects
                            index  = c("property", "year")
    )$coefficients
    
    # run two-way fixed effects with county level data
    coeffmatrix[i,12] <- plm(deforrate6 ~  post*treat, 
                             data   = countypert_df, 
                             method = "within", #fixed effects model
                             effect = "twoway", #county and year fixed effects
                             index  = c("countyid", "year")
    )$coefficients
    
    # run two-way fixed effects with propertylevel data
    coeffmatrix[i,13] <- plm(deforrate7 ~  post*treat, 
                            data   = proppert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #property and year fixed effects
                            index  = c("property", "year")
    )$coefficients
    
    # run two-way fixed effects with county level data
    coeffmatrix[i,14] <- plm(deforrate7 ~  post*treat, 
                             data   = countypert_df, 
                             method = "within", #fixed effects model
                             effect = "twoway", #county and year fixed effects
                             index  = c("countyid", "year")
    )$coefficients
    
    toc()
    rm(year1, countypert_df, proppert_df)
    #end for loop  
  }  
  
  # get distribution information from matrix  
  return(list(colMeans(coeffmatrix), colVars(coeffmatrix)))
  
  assign('coeffmatrix',coeffmatrix, envir=.GlobalEnv)
  
  
  #end function  
}  




