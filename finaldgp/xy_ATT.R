library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
source('defor_DGP.R')

#begin function
xy_ATT <- function(n, nobs, years, min_ATT, max_ATT, base_0, base_1, trend, std_a = 0.1, std_v = 0.25, outcome_var = "y"){
  
  ATT_vals <- seq(from = min_ATT, to = max_ATT, by = .01)
  list <- as.list(ATT_vals)
  
  #preallocate n x ___ matrix
  coeffmatrix <- matrix(nrow = n, ncol = length(list))
  
  b0 = qnorm(base_0, mean = 0, sd = (std_a^2+std_v^2)^.5)
  b1 = qnorm(base_1, mean = 0, sd = (std_a^2+std_v^2)^.5) - b0
  b2 = qnorm(trend + base_0, mean = 0, sd = (std_a^2+std_v^2)^.5) - b0
  
  for(k in list){
    
    b3 = qnorm( pnorm(b0+b1+b2, mean = 0, sd = (std_a^2+std_v^2)^.5) + k , mean = 0, sd = (std_a^2+std_v^2)^.5) - (b0 + b1 + b2)
    
    for(i in 1:n){
      # call defor_dgp function to simulate dataframe, returned as defor_df  
      dgp_results <- defor_DGP(nobs, years, b0, b1, b2, b3, std_a, std_v)
      panels = dgp_results$panels
      ATT = dgp_results$ATT
      DID_estimand = dgp_results$DID_estimand
      panels['outcome'] = panels[outcome_var]
      
      place <- which(list == k)
      
      #mod <- lm(y_it ~  post*treat, data = panels)
      mod <- lm(outcome ~  post*treat, data = panels)
      coeffmatrix[i,place] <- mod$coefficients[4] - ATT

    }
    
    print(k)
    
  }
  
  bias_df <- as.data.frame(cbind(colMeans(coeffmatrix), ATT_vals)) %>%
    mutate(prop = V1/ATT_vals)
  
  
  
  plot <- ggplot(data = bias_df, aes(x = ATT_vals, y = V1)) + 
    geom_line(aes(x = ATT_vals, y = V1, colour = "DID bias"), color = "purple", size =1) +
    #geom_line(aes(x = ATT_vals, y = prop, ymin = 0.001), color = "yellow", size =1) +
    #geom_line(aes(x = ATT_vals, y = prop, ymax = -.001), color = "yellow", size =1) +
    labs(x = "ATT", y = "Bias", caption = paste("Bias as a function of ATT")) + 
    geom_hline(yintercept = 0, linetype = "dashed")
    
  
  outputs = list("bias_df" = bias_df, "plot" = plot)
  
  return(outputs)
  
  #end function  
}  