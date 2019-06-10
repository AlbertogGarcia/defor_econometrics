# function designed to give mean and variance of coefficients when aggregated to county or property level

library(clubSandwich)
library(matrixStats)


#begin function
agg_coeffdist_fcn <- function(n, nobs, years, psize, cellsize){
  
  #preallocate n x 2 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 2)
  
  for(i in 1:n){
    
    
    # call defor_sim function to simulate dataframe, returned as defor_df  
    prop_deforsim(nobs, years, psize, cellsize)  
    
    
    

    
    # run two-way fixed effects with propertylevel data
    coeffmatrix[i,1] <- plm(defor ~  post*treat, 
                            data   = proppert_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #property and year fixed effects
                            index  = c("property", "year")
    )$coefficients
    
    # run two-way fixed effects with county level data
    coeffmatrix[i,2] <- plm(defor ~  post*treat, 
                            data   = countyprop_df, 
                            method = "within", #fixed effects model
                            effect = "twoway", #county and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients
    
    
    
    #end for loop  
  }  
  
  
  # get distribution information from matrix  
  return(list(colMeans(coeffmatrix), colVars(coeffmatrix)))
  
  
  #end function  
}  




