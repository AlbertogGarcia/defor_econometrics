library(ggplot2)
#begin function
DID_bias <- function(n, nobs, years, b0, b1p, b2p, b3){
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  
  for(i in 1:n){
    

   b1n <- -b1p

   b2n <- -b2p
    
     
    # call defor_dgp function to simulate dataframe, returned as defor_df  
    defor_DGP(nobs, years, b0, b1p, b2p, b3)
    
    # run DID dropping deforested pixels
    coeffmatrix[i,1]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
   
    
    defor_DGP(nobs, years, b0, b1p, b2n, b3)
    coeffmatrix[i,2]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
  
    
    defor_DGP(nobs, years, b0, b1n, b2p, b3)
    coeffmatrix[i,3]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
    
    
    defor_DGP(nobs, years, b0, b1n, b2n, b3)
    coeffmatrix[i,4]  <- lm(y_it ~  post*treat, 
                            data = panels
    )$coefficients[4] - ATT
    
    print(i)
    
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
    scale_fill_discrete(breaks=c("++", "+-", "-+", "--"), labels=c("b1+, b2+", "b1+, b2-", "b1-, b2+", "b1-, b2-"))+
    geom_vline(xintercept = 0, linetype = "dashed")+
    #theme(plot.margin = unit(c(1,1,3,1), "cm"))+
    theme(plot.caption = element_text(hjust = 0.5))
    #labs(x= "bias", title = "DID bias with binary outcome", caption = paste("ATT = ",att , "for each of the four density plots. n=", n , "\n",  "The mean and variance are:", "\n", "diff+, trend+:", round(colMeans(coeff_didbias)[1], digits = 4),"var:",round(colVars(coeff_didbias)[1], digits = 6), "\n",  "diff+, trend-:", round(colMeans(coeff_didbias)[2], digits = 4), "var:",round(colVars(coeff_didbias)[2], digits = 6), "\n",  "diff-, trend+:",round(colMeans(coeff_didbias)[3], digits = 4) ,"var:",round(colVars(coeff_didbias)[3], digits = 6),  "\n", "diff-, trend+:",round(colMeans(coeff_didbias)[4], digits = 4) ,"var:",round(colVars(coeff_didbias)[4], digits = 6)))
  pdid
  
  
  #end function  
}  




