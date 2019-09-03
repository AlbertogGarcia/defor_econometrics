library(margins)
#begin function
funcform <- function(n, nobs, years, b0, b1, b2, b3, std_a = 0.1, std_v = 0.25){
  
  #preallocate n x 3 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  
  for(i in 1:n){
    
    
    
    # call defor_sim function to simulate dataframe, returned as defor_df  
    dgp_results <- defor_DGP(nobs, years, b0, b1, b2, b3, std_a, std_v)
    panels = dgp_results$panels
    ATT = dgp_results$ATT
    
    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
    
    
    
    probitcoeff  <- glm(y_it ~  post*treat, 
                        family = binomial(link = "probit"), 
                        data = panels
    )$coefficients
    
    tr <- sum(probitcoeff)
    untr <- probitcoeff[1]+probitcoeff[2]+probitcoeff[3]
    probt <- pnorm(tr)
    probu <- pnorm(untr)
    coeffmatrix[i,2] <- probt-probu - ATT
    
    
    #run logit regression
    
    logitcoeff  <- glm(y_it ~  post*treat, 
                       family = binomial, 
                       data = panels
    )$coefficients
    
    
    
    #need to convert coefficient to probability
    tr <- sum(logitcoeff)
    untr <- logitcoeff[1]+logitcoeff[2]+logitcoeff[3]
    oddst <- exp(tr)
    probt <- oddst / (1 + oddst)
    oddsu <- exp(untr)
    probu <- oddsu / (1 + oddsu)
    coeffmatrix[i,3] <- probt-probu - ATT
    
    
    
    
    
    #end for loop  
  }  
  
  
  
  did_coeff <- as.data.frame(coeffmatrix)
  names(did_coeff)[1] <- paste("DID")
  names(did_coeff)[2] <- paste("Probit")
  names(did_coeff)[3] <- paste("Logit")
  #names(did_coeff)[4] <- paste("Poisson")
  suppressWarnings(did_coeff <- melt(did_coeff, value.name = "bias"))
  
  pdid <- ggplot(data = did_coeff, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("DID", "Probit", "Logit", "Poisson"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "bias", title = "DID bias based on functional form", caption = paste("The mean and variance are:", "\n", "DID:", round(colMeans(coeffmatrix)[1], digits = 4),", var:",round(colVars(coeffmatrix)[1], digits = 6), "\n",  "Probit:", round(colMeans(coeffmatrix)[2], digits = 4), ", var:",round(colVars(coeffmatrix)[2], digits = 6), "\n",  "Logit:",round(colMeans(coeffmatrix)[3], digits = 4) ,", var:",round(colVars(coeffmatrix)[3], digits = 6)))#,  "\n", "Logit:",round(colMeans(coeff_didbias)[4], digits = 4) ,"var:",round(colVars(coeff_didbias)[4], digits = 6)))
  
  outputs = list("plot" = pdid, "did_biases" = did_coeff)
  return(outputs)
  
  #end function  
}  




