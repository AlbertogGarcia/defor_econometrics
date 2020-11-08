# function to see bias from different outcomes as we vary deforestation rates
library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(tictoc)


vary_rate_funcform <- function(n, nobs, years, base_0_start, base_1_start, trend_start, ATT_start, std_a = 0.1, std_v = 0.25){
  
  multipliers <- c(.25, .5, 1, 2.5, 5, 7.5, 10)
  mult_list <- as.list(multipliers)
  
  #preallocate n x ___ matrix
  out_matrix1 <- matrix(nrow = n, ncol = length(mult_list))
  out_matrix2 <- matrix(nrow = n, ncol = length(mult_list))
  out_matrix3 <- matrix(nrow = n, ncol = length(mult_list))
  #out_matrix4 <- matrix(nrow = n, ncol = length(mult_list))
  
  for(k in mult_list){
    tic("loop")
    
    ATT = ATT_start*k
    base_0 = base_0_start*k
    base_1 = base_1_start*k
    trend = trend_start*k
    
    place <- which(mult_list == k)
    
    #compute relevant beta parameters
    std_av = (std_a^2+std_v^2)^.5
    b0 = qnorm(base_0, mean = 0, sd = std_av)
    b1 = qnorm(base_1, mean = 0, sd = std_av) - b0
    b2_0 = qnorm(trend + base_0, mean = 0, sd = std_av) - b0
    b2_1 = qnorm(trend + base_1, mean = 0, sd = std_av) - b0 - b1
    b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_av) + ATT , mean = 0, sd = std_av) - (b0 + b1 + b2_1)
    
    for(i in 1:n){

      # call defor_sim function to simulate dataframe, returned as defor_df  
      dgp_results <- deforestation_DGP(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
      panels = dgp_results$panels
      ATT = dgp_results$ATT
        
      # run DID dropping deforested pixels
      out_matrix1[i,place]  <- (lm(y_it ~  post*treat, 
                                data = panels
      )$coefficients[4] - ATT)/ATT
        
        
        
      probit_me <- probitmfx(formula = y_it ~  post*treat, data = panels)$mfxest[3]
      out_matrix2[i,place] <- (probit_me - ATT)/ATT
        
        
      #run logit regression
        
      logit_me  <- logitmfx(formula = y_it ~  post*treat, data = panels)$mfxest[3]
      out_matrix3[i,place] <- (logit_me - ATT)/ATT
        
        
      # poiss_me <- poissonmfx(formula = y_it ~  post*treat, data = panels)$mfxest[3]
      # out_matrix4[i,place] <- (poiss_me - ATT)/ATT
      #   
        
      #end for loop  
      }
    
    print(k)
    toc()
    
  } # end of list loop
  
  outbias_df1 <- as.data.frame(cbind(colMeans(out_matrix1), multipliers)) 
  outbias_df2 <- as.data.frame(cbind(colMeans(out_matrix2), multipliers)) 
  outbias_df3 <- as.data.frame(cbind(colMeans(out_matrix3), multipliers)) 
  
  colors <- c("Linear Probability Model" = "red", 
              "Probit" = "yellow",
              "Logit" = "blue")
  
  plot <- ggplot() + 
    geom_line(data = outbias_df1, aes(x = multipliers, y = V1, color="Linear Probability Model"), size =1) +
    geom_line(data = outbias_df2, aes(x = multipliers, y = V1, color="Probit"), size =1) +
    geom_line(data = outbias_df3, aes(x = multipliers, y = V1, color="Logit"), size =1) +
    #geom_line(data = outbias_df4, aes(x = multipliers, y = V1, color="Poisson"), size =1) +
    labs(x = "Ratio of landscape parameters relative to guiding example parameterization", y = "Bias as a proportion of ATT") +
    scale_color_manual(values = colors)+
    geom_hline(yintercept = 0, linetype = "dashed")
  
  outputs = list("outbias_df1" = outbias_df1, "outbias_df2" = outbias_df2, "outbias_df3" = outbias_df3,"plot" = plot)
  
  
} # end of function
