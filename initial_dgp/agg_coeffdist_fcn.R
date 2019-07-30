# function designed to give mean and variance of coefficients when aggregated to county or property level

library(clubSandwich)
library(matrixStats)
library(DataCombine)
library(plm)
library(IDPmisc)
library(tictoc)
#begin function
agg_coeffdist_fcn <- function(n, nobs, years, psize, cellsize, std_p, att){
  
  #preallocate n x 2 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  areamatrix <- matrix(nrow = n, ncol = 2)
  for(i in 1:n){
    
    tic("loop")
    # call defor_sim function to simulate dataframe, returned as defor_df  
    prop_deforsim(nobs, years, psize, cellsize, std_p, att)  
    
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
    
    ##### constructing outcome variables #####
    proppert_df$forcovariate <- log(proppert_df$forshare * proppert_df$parea)
    
    proppert_df$deforrate1 <- ((proppert_df$forsharelag - proppert_df$forshare) / proppert_df$forsharelag)
    
    
    countypert_df$deforrate1 <- ((countypert_df$forsharelag- countypert_df$forshare) / countypert_df$forsharelag)
    
    proppert_df$deforrate2 <- log(1+ proppert_df$defor*proppert_df$parea)
    #remove infinite values

    proppert_df <- subset(proppert_df, select = -c(geometry))
    proppert_df <- 
    proppert_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    countypert_df <- subset(countypert_df, select = -c(geometry))
    countypert_df <- 
      countypert_df %>% 
      filter_all(all_vars(!is.infinite(.)))
    
    
    # run two-way fixed effects with county level data
    coeffmatrix[i,1] <- plm(deforrate1 ~  post*treat, 
                            data   = countypert_df, 
                            weights=(countypert_df$carea),
                            method = "within", #fixed effects model
                            effect = "twoway", #county and year fixed effects
                            index  = c("countyid", "year")
    )$coefficients - att
    
    # run two-way fixed effects with propertylevel data
    coeffmatrix[i,2] <- plm(deforrate1 ~  post*treat, 
                            data   = proppert_df, 
                            weights=(proppert_df$parea),
                            method = "within", #fixed effects model
                            effect = "twoway", #property and year fixed effects
                            index  = c("property", "year")
    )$coefficients - att
    
    
    coeffmatrix[i,3] <- plm(deforrate2 ~  post*treat + forcovariate, 
                            data   = proppert_df, 
                            weights=(proppert_df$parea),
                            method = "within", #fixed effects model
                            effect = "twoway", #property and year fixed effects
                            index  = c("property", "year")
    )$coefficients - att
    
    areamatrix[i,1] <- mean(countypert_df$carea)
    areamatrix[i,2]<- mean(proppert_df$parea)
    
    toc()
    rm(countypert_df, proppert_df)
    #end for loop  
  }  
  
  # get distribution information from matrix  
  assign('coeff_agg',coeffmatrix, envir=.GlobalEnv)
  assign('areamatrix',areamatrix, envir=.GlobalEnv)
  return(list(colMeans(coeffmatrix), colVars(coeffmatrix)))
  
  
  
  
  #end function  
}  




