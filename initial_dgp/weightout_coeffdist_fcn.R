library(clubSandwich)
library(matrixStats)
library(DataCombine)
library(plm)
library(IDPmisc)
library(tictoc)
#begin function
weightout_coeffdist_fcn <- function(n, nobs, years, psize, cellsize){
  
  #preallocate n x _ matrix
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  
  for(i in 1:n){
    
    tic("loop")
    # call defor_simkeep function to simulate dataframe, returned as deforkeep_df
    
    
    
    defor_sim(nobs, years)
    property_fcn(defor_df, psize,cellsize)
    
    
    # lagging deforestation and forest share variable
    
    propertylevel_df <- propertylevel_df[order(propertylevel_df$property, propertylevel_df$year),]
    propertylevel_df <- slide(propertylevel_df, Var = "defor", GroupVar = "property", NewVar = "deforlag",
                          slideBy = -1, reminder = FALSE)
    
    
    ##### creating forested share variable #####
    propertylevel_df$forshare <- 1 -  propertylevel_df$defor
    
    propertylevel_df$forsharelag <- 1 -  propertylevel_df$deforlag
    
    
    # #### create baseline defor and forshare vars 
    # #county level
    # year1 <- subset(gridlevel_df, year== 1)
    # colnames(year1)[colnames(year1)=="defor"] <- "defor0"
    # year1 <- subset(year1, select = c(countyid, defor0))
    # #st_geometry(year1) <- NULL
    # gridlevel_df <-  merge(gridlevel_df, year1, by = "countyid")
    # 
    # gridlevel_df$forshare0 <- 1 -  gridlevel_df$defor0
    
    
    
    ##### constructing various outcome variables #####
    propertylevel_df$deforrate1 <- (propertylevel_df$forshare - propertylevel_df$forsharelag)/ propertylevel_df$parea
    
    propertylevel_df$deforrate2 <- (propertylevel_df$parea*(propertylevel_df$defor - propertylevel_df$deforlag)) / propertylevel_df$forsharelag
    
    propertylevel_df$deforrate3 <- ((propertylevel_df$forsharelag- propertylevel_df$forshare) / propertylevel_df$forsharelag)
    
    
    #remove any infinite values
      propertylevel_df <- subset(propertylevel_df, select = -c(geometry))
    propertylevel_df <- 
      propertylevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    
    
    # run two-way fixed effects with outcome 1 
    coeffmatrix[i,1] <- plm(deforrate1 ~  post*treat, 
                            data   = propertylevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("property", "year")
    )$coefficients
    
    coeffmatrix[i,2] <- plm(deforrate2 ~  post*treat, 
                            data   = propertylevel_df, 
                            weights=(propertylevel_df$parea),
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("property", "year")
                            
    )$coefficients
    
    
    coeffmatrix[i,3] <- plm(deforrate2 ~  post*treat, 
                            data   = propertylevel_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #countyid and year fixed effects
                            index  = c("property", "year")
    )$coefficients
    
   
    rm(propertylevel_df)
    toc()
    #end for loop
  }  
  
  # get distribution information from matrix  
  assign('coeff_weight',coeffmatrix, envir=.GlobalEnv)
  return(list(colMeans(coeffmatrix), colVars(coeffmatrix)))
  
  
  #end function  
}  




