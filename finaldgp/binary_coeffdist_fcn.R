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
    ATT = dgp_results$ATT
    DID_estimand = dgp_results$DID_estimand
    
    DIFF = pnorm(b0+b1, 0, (std_a^2+std_v^2)^.5) - pnorm(b0, 0, (std_a^2+std_v^2)^.5)

    
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - DID_estimand
    
    # run two-way fixed effects    
    coeffmatrix[i, 2] <- plm(y_it ~  post*treat, 
                             data   = panels, 
                             method = "within", #fixed effects model
                             effect = "twoway", #unit and year fixed effects
                             index  = c("pixels", "year")
    )$coefficients - DID_estimand
    
    
    
    
    # DID keeping variables
    
    coeffmatrix[i, 3] <- lm(defor ~  post*treat, 
                           data = panels
    )$coefficients[4] - DID_estimand
    
    # run two-way fixed effects    
    coeffmatrix[i, 4] <- plm(defor ~  post*treat, 
                       data   = panels, 
                       method = "within", #fixed effects model
                       effect = "twoway", #unit and year fixed effects
                       index  = c("pixels", "year")
    )$coefficients - DID_estimand
    
    
  #end for loop  
  }  
  #assign('coeffmatrix',coeffmatrix, envir=.GlobalEnv)
  
# get distribution information from matrix  

  b_coeff <- as.data.frame(coeffmatrix)
  actual <- rep(DID_estimand, times = n)
  
  names(b_coeff)[1] <- paste("DID1")
  names(b_coeff)[2] <- paste("tw1")
  
  names(b_coeff)[3] <- paste("DID2")
  names(b_coeff)[4] <- paste("tw2")
  suppressWarnings(cbias <- melt(b_coeff, value.name = "bias"))
  
  plot = ggplot(data = cbias, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("DID1", "tw1", "DID2", "tw2"), labels=c("DID dropping obs", "2way FE dropping obs", "DID keeping obs", "2way FE keeping obs"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #geom_vline(aes(xintercept= (DID_estimand - ATT), color="DID estimand - ATT"), linetype="dashed", color = "green")+
    #geom_vline(aes(xintercept= (DIFF), color="2way FE bias from group DIFF"), linetype="dashed", color = "red")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", title = "Bias with binary outcome", caption = paste("The mean and RMSE are:", "\n", "DID dropping obs:", round(colMeans(coeffmatrix)[1], digits = 4),"RMSE:",round(rmse(actual, coeffmatrix[1]), digits = 5), "\n", 
                                                                        "2way FE dropping obs:", round(colMeans(coeffmatrix)[2], digits = 4), "RMSE:",round(rmse(actual, coeffmatrix[2]), digits = 5), "\n",  
                                                                        "DID keeping obs:",round(colMeans(coeffmatrix)[3], digits = 4) ,"RMSE:",round(rmse(actual, coeffmatrix[3]), digits = 5),  "\n", 
                                                                        "2way FE keeping obs:",round(colMeans(coeffmatrix)[4], digits = 4) ,"RMSE:",round(rmse(actual, coeffmatrix[4]), digits = 5))
    )
  
  model <- c("DID_drop", "TWFE_drop", "DID_keep", "TWFE_keep")
  bias <- c(colMeans(coeffmatrix)[1], colMeans(coeffmatrix)[2], colMeans(coeffmatrix)[3], colMeans(coeffmatrix)[4])
  RMSE <- c(rmse(actual, coeffmatrix[1]), rmse(actual, coeffmatrix[2]), rmse(actual, coeffmatrix[3]), rmse(actual, coeffmatrix[4]))
  binary_results <- data.frame(model, bias, RMSE)
  outputs = list("plot" = plot, "did_biases" = b_coeff, "binary_results" = binary_results)
  return(outputs)
  
#end function  
}  

  
  
  
  