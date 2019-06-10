# this function's intention is to provide distributional parameters for the various specifications that have the typical binary outcome variable
library(clubSandwich)
library(matrixStats)


#begin function
binary_coeffdist_fcn <- function(n, nobs, years){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  
  for(i in 1:n){
  
  
    # call defor_sim function to simulate dataframe, returned as defor_df  
    defor_sim(nobs, years)  

   
    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(defor ~  post*treat, 
                   data = defor_df
    )$coefficients[4]
    
    # run two-way fixed effects    
    coeffmatrix[i,2] <- plm(defor ~  post*treat, 
                     data   = defor_df, 
                     method = "within", #fixed effects model
                     effect = "twoway", #unit and year fixed effects
                     index  = c("idx", "year")
    )$coefficients
    
    # turn dropped outcome observations to 1
    defor_df$defor[is.na(defor_df$defor)] <- 1
    
    # DID keeping variables
    
    coeffmatrix[i,3] <- lm(defor ~  post*treat, 
                           data = defor_df
    )$coefficients[4]
    
    # run two-way fixed effects    
    coeffmatrix[i,4] <- plm(defor ~  post*treat, 
                       data   = defor_df, 
                       method = "within", #fixed effects model
                       effect = "twoway", #unit and year fixed effects
                       index  = c("idx", "year")
    )$coefficients
    
    rm(defor_df)
    
  #end for loop  
  }  
  
  
# get distribution information from matrix  
  return(list(colMeans(coeffmatrix), colVars(coeffmatrix)))

  
#end function  
}  
  
  
  
  
  