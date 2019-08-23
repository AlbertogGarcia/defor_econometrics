# coverage function 
library(plm)
library(clubSandwich)

coverage <- function(n, nobs, years, b0, b1, b2, b3){
  
  covmat <- matrix(nrow = n, ncol = 3)
  
  for(i in 1:n){
    
    defor_DGP(nobs, years, b0, b1, b2, b3)
    
    DID <- lm(y_it ~  post*treat, 
       data = panels)
    
    
    vcov <- vcovHC(DID, type = "HC0")
    cluster_se    <- sqrt(diag(vcov))
    
    covmat[i,1] <- ATT
    
    covmat[i,2] <- DID$coefficients[4] - 1.96 * cluster_se[4]
    
    covmat[i,3] <- DID$coefficients[4] + 1.96*  cluster_se[4]
    
      
  }
  
  
  coverage_df <- as.data.frame(covmat)
  
  coverage_df$indic <- ifelse(coverage_df$V1 < coverage_df$V2 | coverage_df$V1 > coverage_df$V3 , 0, 1)
  
  print(mean(coverage_df$indic))
}