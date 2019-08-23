library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
#begin function
quickmonteyit <- function(n, nobs, years, b0, b1, b2, b3){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 1)
  
  for(i in 1:n){
    
    
    
    # call defor_dgp function to simulate dataframe, returned as defor_df  
    defor_DGP(nobs, years, b0, b1, b2, b3)
    
    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
    
    print(i)
    
    #end for loop  
  }  
  
  #assign('coeff_didbias',coeffmatrix, envir=.GlobalEnv)
  
  did_coeff <- as.data.frame(coeffmatrix)
  actual <- rep(ATT, times = n)
  
  ggplot() +
    geom_density(data = did_coeff , aes(x = V1), alpha = .2, fill='red')+
    geom_vline(data = did_coeff, xintercept = mean(did_coeff$V1), color = 'red')+
    geom_vline(data = did_coeff, xintercept = 0, linetype = "dashed")+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "Bias", caption = paste("Mean:", round(mean(did_coeff$V1), digits = 4),
                                    ", RMSE:", round(rmse(actual, did_coeff$V1), digits = 4)) 
    )
  
  
  
  #end function  
}  
