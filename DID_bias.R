
#begin function
DID_bias <- function(n, nobs, years, att){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  
  for(i in 1:n){
    
    diff1 <- .2
    trend1 <- .1
    diff2 <- -diff1
    trend2 <- -trend1  
    # call defor_sim function to simulate dataframe, returned as defor_df  
    defor_sim(nobs, years, att, diff1, trend1)
    
    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(defor ~  post*treat, 
                            data = defor_df
    )$coefficients[4] - att
   
    
    defor_sim(nobs, years, att, diff1, trend2)
    coeffmatrix[i,2]  <- lm(defor ~  post*treat, 
                            data = defor_df
    )$coefficients[4] - att
  
    
    defor_sim(nobs, years, att, diff2, trend1)
    coeffmatrix[i,3]  <- lm(defor ~  post*treat, 
                            data = defor_df
    )$coefficients[4] - att
    
    
    defor_sim(nobs, years, att, diff2, trend2)
    coeffmatrix[i,4]  <- lm(defor ~  post*treat, 
                            data = defor_df
    )$coefficients[4] - att
    

    
    #end for loop  
  }  
  assign('coeff_didbias',coeffmatrix, envir=.GlobalEnv)
  
  did_coeff <- as.data.frame(coeff_didbias)
  names(did_coeff)[1] <- paste("++")
  names(did_coeff)[2] <- paste("+-")
  names(did_coeff)[3] <- paste("-+")
  names(did_coeff)[4] <- paste("--")
  suppressWarnings(did_coeff <- melt(did_coeff, value.name = "bias"))
  
  pdid <- ggplot(data = did_coeff, aes(x = bias, fill=variable)) +
    geom_density(alpha = .2) +
    guides(fill=guide_legend(title=NULL))+
    scale_fill_discrete(breaks=c("++", "+-", "-+", "--"), labels=c("diff+, trend+", "diff+, trend-", "diff-, trend+", "diff-, trend-"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))+
    labs(x= "bias", title = "DID bias with binary outcome", caption = paste("ATT = ",att , "for each of the four density plots. n=", n , "\n",  "The mean and variance are:", "\n", "diff+, trend+:", round(colMeans(coeff_didbias)[1], digits = 4),"var:",round(colVars(coeff_didbias)[1], digits = 6), "\n",  "diff+, trend-:", round(colMeans(coeff_didbias)[2], digits = 4), "var:",round(colVars(coeff_didbias)[2], digits = 6), "\n",  "diff-, trend+:",round(colMeans(coeff_didbias)[3], digits = 4) ,"var:",round(colVars(coeff_didbias)[3], digits = 6),  "\n", "diff-, trend+:",round(colMeans(coeff_didbias)[4], digits = 4) ,"var:",round(colVars(coeff_didbias)[4], digits = 6)))
  pdid
  
  
  #end function  
}  




