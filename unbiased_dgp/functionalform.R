library(margins)
library(mfx)
#begin function
functionalform <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0, std_v = 0.25){
  
  #preallocate n x 3 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  
  for(i in 1:n){
    
    
    
    # call defor_sim function to simulate dataframe, returned as defor_df  
    dgp_results <- deforestation_DGP(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
    panels = dgp_results$panels
    ATT = dgp_results$ATT
    
    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
    
    
    
    probit_me <- probitmfx(formula = y_it ~  post*treat, data = panels)$mfxest[3]
    coeffmatrix[i,2] <- probit_me - ATT
    
    
    #run logit regression
    
    logit_me  <- logitmfx(formula = y_it ~  post*treat, data = panels)$mfxest[3]
    coeffmatrix[i,3] <- logit_me - ATT
    
    
    # poiss_me <- poissonmfx(formula = y_it ~  post*treat, data = panels)$mfxest[3]
    # coeffmatrix[i,4] <- poiss_me - ATT
    
    
    #end for loop  
  }  
  
  actual <- rep(0, times = n)
  
  did_coeff <- as.data.frame(coeffmatrix)
  names(did_coeff)[1] <- paste("LPM")
  names(did_coeff)[2] <- paste("Probit")
  names(did_coeff)[3] <- paste("Logit")
  names(did_coeff)[4] <- paste("Poisson")
  suppressWarnings(did_coeff <- melt(did_coeff, value.name = "bias"))
  
  pdid <- ggplot(data = did_coeff, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("LPM", "Probit", "Logit"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "bias", title = "DID bias based on functional form", 
         caption = paste("The mean and variance are:", "\n", 
                         "LPM:"    , round(colMeans(coeffmatrix)[1], digits = 4), "RMSE:", round(rmse(actual, coeffmatrix[1]), digits = 5), "\n",  
                         "Probit:" , round(colMeans(coeffmatrix)[2], digits = 4), "RMSE:", round(rmse(actual, coeffmatrix[2]), digits = 5), "\n",  
                         "Logit:"  , round(colMeans(coeffmatrix)[3], digits = 4), "RMSE:", round(rmse(actual, coeffmatrix[3]), digits = 5)
                         #, "\n", "Poisson:", round(colMeans(coeffmatrix)[4], digits = 4), "RMSE:", round(rmse(actual, coeffmatrix[4]), digits = 5)
         )
    )
  
  
  outputs = list("plot" = pdid, "did_biases" = did_coeff)
  return(outputs)
  
  #end function  
}  




