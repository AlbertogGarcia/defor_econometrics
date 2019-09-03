# this function's intention is to provide distributional parameters for the various specifications that have the typical binary outcome variable
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
source('defor_DGP.R')
#begin function

binary_coeffdist_fcn <- function(n, nobs, years, b0, b1, b2, b3, std_a = 0.1, std_v = 0.25){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  
  for(i in 1:n){
    
  
  
    # call defor_sim function to simulate dataframe, returned as panels  
    dgp_results <- defor_DGP(nobs, years, b0, b1, b2, b3, std_a, std_v)
    panels = dgp_results$panels

    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                   data = panels
    )$coefficients[4] - ATT
    
    # run two-way fixed effects    
    coeffmatrix[i, 2] <- plm(y_it ~  post*treat, 
                     data   = panels, 
                     method = "within", #fixed effects model
                     effect = "twoway", #unit and year fixed effects
                     index  = c("pixels", "year")
    )$coefficients - ATT
    
    
    # DID keeping variables
    
    coeffmatrix[i, 3] <- lm(defor ~  post*treat, 
                           data = panels
    )$coefficients[4] - ATT
    
    # run two-way fixed effects    
    coeffmatrix[i, 4] <- plm(defor ~  post*treat, 
                       data   = panels, 
                       method = "within", #fixed effects model
                       effect = "twoway", #unit and year fixed effects
                       index  = c("pixels", "year")
    )$coefficients - ATT
    
    
  #end for loop  
  }  
  assign('coeff_bin',coeffmatrix, envir=.GlobalEnv)
  
# get distribution information from matrix  

  b_coeff <- as.data.frame(coeff_bin)
  
  actual <- rep(ATT, times = n)
  
  p1 <- ggplot() +
    geom_density(data = b_coeff , aes(x = V1), alpha = .2, fill="#29CD44") +
    geom_density(data = b_coeff , aes(x = V2), alpha = .2, fill="#FF6655") +
    geom_density(data = b_coeff , aes(x = V3), alpha = .2, fill="#0000CC") +
    geom_density(data = b_coeff , aes(x = V4), alpha = .2, fill="#0000CC") +
    geom_vline(data = b_coeff, xintercept = mean(b_coeff$V1),linetype = "dashed")+
    geom_vline(data = b_coeff, xintercept = mean(b_coeff$V2),linetype = "dashed")+
    geom_vline(data = b_coeff, xintercept = mean(b_coeff$V3),linetype = "dashed")+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", title = "Bias with binary outcome", caption = paste("The mean and RMSE are:", "\n", "DID dropping obs:", round(colMeans(coeff_bin)[1], digits = 4),"RMSE:",round(rmse(actual, coeff_bin[1]), digits = 5), "\n", 
                                                                        "2way FE dropping obs:", round(colMeans(coeff_bin)[2], digits = 4), "RMSE:",round(rmse(actual, coeff_bin[2]), digits = 5), "\n",  
                                                                        "DID keeping obs:",round(colMeans(coeff_bin)[3], digits = 4) ,"RMSE:",round(rmse(actual, coeff_bin[3]), digits = 5),  "\n", 
                                                                        "2way FE keeping obs:",round(colMeans(coeff_bin)[4], digits = 4) ,"RMSE:",round(rmse(actual, coeff_bin[4]), digits = 5)))
  p1
#end function  
}  

  
  
  
  
