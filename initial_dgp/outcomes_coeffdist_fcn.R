library(clubSandwich)
library(matrixStats)
library(DataCombine)
library(plm)
library(IDPmisc)
library(tictoc)
#begin function
outcomes_coeffdist_fcn <- function(n, nobs, years, cellsize){
  
  #preallocate n x _ matrix
  coeffmatrix <- matrix(nrow = n, ncol = 5)
  
  for(i in 1:n){
    
    tic("loop")
    # call defor_simkeep function to simulate dataframe, returned as deforkeep_df
    
    
    
    defor_sim(nobs, years)
    grid_fcn(defor_df, cellsize)
    #prop_deforsim(10000,3,1200,25)  
    ############### gather area for each property and county ##################
    
    
    
    
    # lagging deforestation and forest share variable
    
    gridlevel_df <- gridlevel_df[order(gridlevel_df$countyid, gridlevel_df$year),]
    gridlevel_df <- slide(gridlevel_df, Var = "defor", GroupVar = "countyid", NewVar = "deforlag",
                           slideBy = -1, reminder = FALSE)
   
    
    ##### creating forested share variable #####
    gridlevel_df$forshare <- 1 -  gridlevel_df$defor
  
    gridlevel_df$forsharelag <- 1 -  gridlevel_df$deforlag

    
    #### create baseline defor and forshare vars 
    #county level
    year1 <- subset(gridlevel_df, year== 1)
    colnames(year1)[colnames(year1)=="defor"] <- "defor0"
    year1 <- subset(year1, select = c(countyid, defor0))
    #st_geometry(year1) <- NULL
    gridlevel_df <-  merge(gridlevel_df, year1, by = "countyid")
   
    gridlevel_df$forshare0 <- 1 -  gridlevel_df$defor0

    
    
    ##### constructing various outcome variables #####
    gridlevel_df$deforrate1 <- gridlevel_df$forshare / gridlevel_df$forsharelag
    
    gridlevel_df$deforrate2 <- ((gridlevel_df$forsharelag- gridlevel_df$forshare) / gridlevel_df$forsharelag)
    
    gridlevel_df$deforrate3 <- gridlevel_df$forshare / gridlevel_df$forshare0

    gridlevel_df$deforrate4 <- gridlevel_df$defor / gridlevel_df$forshare0

    gridlevel_df$deforrate5 <- gridlevel_df$defor / gridlevel_df$forsharelag
    
    #remove any infinite values
    gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
    gridlevel_df <- 
      gridlevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    
    
    # run two-way fixed effects with outcome 1 
    coeffmatrix[i,1] <- plm(deforrate1 ~  post*treat, 
                            data   = gridlevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    coeffmatrix[i,2] <- plm(deforrate2 ~  post*treat, 
                            data   = gridlevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    
    coeffmatrix[i,3] <- plm(deforrate3 ~  post*treat, 
                            data   = gridlevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    coeffmatrix[i,4] <- plm(deforrate4 ~  post*treat, 
                            data   = gridlevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    coeffmatrix[i,5] <- plm(deforrate5 ~  post*treat,
                            data   = gridlevel_df,
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    
    
    rm(year1, gridlevel_df)
    toc()
    #end for loop
  }  
  
  # get distribution information from matrix  
  assign('coeff_out',coeffmatrix, envir=.GlobalEnv)
  return(list(colMeans(coeffmatrix), colVars(coeffmatrix)))
  
  
  #end function  
}  




