library(margins)
#begin function
functionalforms <- function(n, nobs, years, att){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 3)
  
  for(i in 1:n){
    
    diff <- .2
    trend <-  -.1
    
    
    # call defor_sim function to simulate dataframe, returned as defor_df  
    defor_sim(nobs, years, att, diff, trend)
    
    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(defor ~  post*treat, 
                            data = defor_df
    )$coefficients[4] - att
    
    
    
    probitcoeff  <- glm(defor ~  post*treat, 
                            family = binomial(link = "probit"), 
                            data = defor_df
    )$coefficients
    
    tr <- sum(probitcoeff)
    untr <- probitcoeff[1]+probitcoeff[2]+probitcoeff[3]
    probt <- pnorm(tr)
    probu <- pnorm(untr)
    coeffmatrix[i,2] <- probt-probu - att
    
    
    #run logit regression
   
    logitcoeff  <- glm(defor ~  post*treat, 
                             family = binomial, 
                             data = defor_df
    )$coefficients
    
    
    
    #need to convert coefficient to probability
    tr <- sum(logitcoeff)
    untr <- logitcoeff[1]+logitcoeff[2]+logitcoeff[3]
    oddst <- exp(tr)
    probt <- oddst / (1 + oddst)
    oddsu <- exp(untr)
    probu <- oddsu / (1 + oddsu)
    coeffmatrix[i,3] <- probt-probu - att
    
    
    

    
    #end for loop  
  }  
  assign('coeff_didbias',coeffmatrix, envir=.GlobalEnv)
  
  did_coeff <- as.data.frame(coeff_didbias)
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
    labs(x= "bias", title = "DID bias based on functional form", caption = paste("ATT = ",att , "for each of the four density plots. n=", n , "\n",  "The mean and variance are:", "\n", "DID:", round(colMeans(coeff_didbias)[1], digits = 4),", var:",round(colVars(coeff_didbias)[1], digits = 6), "\n",  "Probit:", round(colMeans(coeff_didbias)[2], digits = 4), ", var:",round(colVars(coeff_didbias)[2], digits = 6), "\n",  "Logit:",round(colMeans(coeff_didbias)[3], digits = 4) ,", var:",round(colVars(coeff_didbias)[3], digits = 6)))#,  "\n", "Logit:",round(colMeans(coeff_didbias)[4], digits = 4) ,"var:",round(colVars(coeff_didbias)[4], digits = 6)))
  pdid
  
  
  #end function  
}  




