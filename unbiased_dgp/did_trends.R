library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggpubr)
library(plm)
library(Metrics)
source('DGP_counterfactual.R')

#begin function
did_trends <- function(nobs, years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25){
  
  # call defor_dgp function to simulate dataframe, returned as defor_df  
  dgp_results <- DGP_counterfactual(nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
  panels = dgp_results$panels
  panels_cf = dgp_results$panels_cf %>%
    filter(treat ==1)
  ATT = dgp_results$ATT
  
  panels_t <- panels %>%
    filter(treat == 1)
  panels_u <- panels %>%
    filter(treat == 0)
  
  colors <- c("treated, observed" = "#a1d76a", 
              "untreated, observed" = "light blue",
              "treated, counterfactual" = "black",
              "start of intervention" = "black")
  
  plot <- ggplot()+
    stat_summary(data = panels_t, aes(x = year, y = y, color = "treated, observed"), geom = 'line', size =1.5)+
    stat_summary(data = panels_u, aes(x = year, y = y, color = "untreated, observed"), geom = 'line', size =1.5)+
    stat_summary(data = panels_cf, aes(x = year, y = y_cf, color = "treated, counterfactual"), geom = 'line', size =1.5, linetype = "dashed")+
    geom_vline(xintercept = years, aes(color = "start of intervention"), linetype = "dashed")+
    scale_color_manual(values = colors)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "gray95"), 
          legend.title = element_blank()
          )
  
  outputs = list("plot" = plot)
  return(outputs)
  
  #end function  
}  




