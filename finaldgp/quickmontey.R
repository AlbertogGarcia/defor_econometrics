library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
source('defor_DGP.R')

#begin function
quickmontey2 <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25, outcome_var = "y"){
  
  #preallocate n x 1 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 1)

  for(i in 1:n){
    # call defor_dgp function to simulate dataframe, returned as defor_df  
    dgp_results <- defor_DGP2(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
    panels = dgp_results$panels
    ATT = dgp_results$ATT
    #DID_estimand = dgp_results$DID_estimand
    panels['outcome'] = panels[outcome_var]
    
    
    # run DID dropping deforested pixels
    # mod  <- plm(outcome ~ post*treat, effect = "twoways",
    #             index = c("pixels", "year"), data = panels)
    # coeffmatrix[i,1]  <- mod$coefficients[1] - ATT
    
    mod <- lm(outcome ~  post*treat, data = panels)
    coeffmatrix[i,1] <- mod$coefficients[4] - ATT
    print(i)
    
    #end for loop  
  }  
  
  #assign('coeff_didbias',coeffmatrix, envir=.GlobalEnv)
  
  did_coeff <- as.data.frame(coeffmatrix)
  actual <- rep(ATT, times = n)
  
  plot = ggplot() +
    geom_density(data = did_coeff , aes(x = V1), alpha = .2, fill="#29CD44")+
    #geom_vline(aes(xintercept= (DID_estimand - ATT), color="DID estimand - ATT"), linetype="dashed")+
    geom_vline(data = did_coeff, xintercept = 0, linetype = "dashed")+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Mean:", round(mean(did_coeff$V1), digits = 4),
                                    ", RMSE:", round(rmse(actual, did_coeff$V1), digits = 4)) 
    )
  
  
  
  outputs = list("plot" = plot, "did_biases" = did_coeff)
  return(outputs)
  
  #end function  
}  




