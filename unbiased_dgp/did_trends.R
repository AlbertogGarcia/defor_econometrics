library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggpubr)
library(plm)
library(Metrics)
library(grid)
library(pBrackets)
source('DGP_counterfactual.R')

#begin function
did_trends <- function(years, b0, b1, b2_0, b2_1, b3, std_a = 0.1, std_v = 0.25){
  
  # call defor_dgp function to simulate dataframe, returned as defor_df  
  dgp_results <- DGP_counterfactual(1000000, years, b0, b1, b2_0, b2_1, b3, std_a, std_v)
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
  
  bracketsGrob <- function(...){
    l <- list(...)
    e <- new.env()
    e$l <- l
    grid:::recordGrob(  {
      do.call(grid.brackets, l)
    }, e)
  }
  
  
  b1 <- bracketsGrob(.08, 0.19, .08, .92, h=-0.02, lwd=1, col="black")
  b2 <- bracketsGrob(.9, 0.07, .9, .54, h=0.02, lwd=1, col="black")
  b3 <- bracketsGrob(.63, 0.57, .63, .81, h=-0.02, lwd=1, col="red")
  
  
  plot <- ggplot()+
    stat_summary(data = panels_t, aes(x = year, y = y, color = "treated, observed"), geom = 'line', size =1.5)+
    stat_summary(data = panels_u, aes(x = year, y = y, color = "untreated, observed"), geom = 'line', size =1.5)+
    stat_summary(data = panels_cf, aes(x = year, y = y_cf, color = "treated, counterfactual"), geom = 'line', size =1.5, linetype = "dashed")+
    geom_vline(xintercept = years, aes(color = "start of intervention"), linetype = "dashed")+
    scale_color_manual(values = colors)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "gray95"),
          legend.title = element_blank(),
          #plot.margin = unit(c(10,10,10,10), "lines")
          )+
    # scale_y_continuous("",breaks=c(.025,.04),labels=c("High types", "test"), position = "right" )+
    # scale_y_continuous("",breaks=c(.035),labels=c("Low types"))+
    annotation_custom(b1)+
    annotation_custom(b2)+
    annotation_custom(b3)+
    annotation_custom(grob = textGrob(label = "pre-treatment difference", x = 0.115, y = .555, hjust = 0))+
    annotation_custom(grob = textGrob(label = "post-treatment difference", x = .54, y = .305, hjust = 0))+
    annotation_custom(grob = textGrob(label = "estimated ATT", x = .68, y = .69, hjust = 0))+
    ylab(label = "deforestation rate")
    
  plot
  
  outputs = list("plot" = plot)
  return(outputs)
  
  #end function  
}  




