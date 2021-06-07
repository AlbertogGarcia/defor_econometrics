#this function generates landscape maps using the same dgp as aggregation_method fcn

library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(DeclareDesign)

select <- dplyr::select

source(here::here("unbiased_dgp", 'full_landscape.R'))

#begin function
landscape_map <- function(nobs, years, b0, b1, b2, b3, std_a = 0.1, std_v = 0.25, std_p = .1, cellsize, ppoints, cpoints){
  
  countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  control_area = countyscape$control_area
  intervention_area = countyscape$intervention_area
  p_bounds = countyscape$p_bounds
  c_bounds = countyscape$c_bounds
  
  
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  pixloc <- pixloc_df
  
  Nobs <- length(pixloc$treat)  
  panels <- fabricate(
    pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
    year = add_level(N = (years*2), nest = FALSE),
    obs = cross_levels(
      by = join(pixels, year),
      post = ifelse(year > years, 1, 0),
      v_it = rnorm(N, 0, std_v),
      ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
      ystar_counterfactual = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
    )
  )
  
  
  #generate random 
  error_table <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, sd = std_p))
  
  panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
  
  panels <- panels %>%
    inner_join(pixloc, by = c("pixels", "treat")) %>%
    inner_join(error_table, by = "property") %>%
    mutate(ystar = ystar + p_err, 
           y = (ystar > 0)*1, 
           ystar_counterfactual = ystar_counterfactual+p_err,
           y_counterfactual = (ystar_counterfactual > 0)*1) 
  
  panels_counterfactual <- panels
  
  #need to determine which year deforestation occurred
  year_df <- panels %>%
    dplyr::select(pixels, year, y) %>%
    dcast(pixels ~ year , value.var = "y")
  
  rownames(year_df) <- year_df$pixels
  
  year_df <- year_df %>%
    dplyr::select(- pixels)
  
  #creating variable for the year a pixel is deforested
  not_defor <- rowSums(year_df)<1 *1
  defor_year <- max.col(year_df, ties.method = "first") 
  defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
  defor_df <- tibble::rownames_to_column(defor_df)
  names(defor_df)[1] <- paste("pixels")
  
  panels <- defor_df %>%
    dplyr::select(pixels, defor_year) %>%
    inner_join(panels, by = "pixels")
  
  cols.num <- c("pixels", "grid", "property", "county", "year")
  panels[cols.num] <- sapply(panels[cols.num],as.numeric)
  
  panels <- panels %>%
    mutate(indic = year - defor_year,
           defor = ifelse(indic > 0, 1, y),
           y_it = ifelse(indic > 0, NA, y) )
           
  
  ##################################################################
  
  ################   creating counterfactual landscape     #########################
  
  ######################################################################
  year_counterfactual <- panels_counterfactual %>%
    dplyr::select(pixels, year, y_counterfactual) %>%
    dcast(pixels ~ year , value.var = "y_counterfactual")
  
  rownames(year_counterfactual) <- year_counterfactual$pixels
  
  year_counterfactual <- year_counterfactual %>%
    dplyr::select(- pixels)
  
  #creating variable for the year a pixel is deforested
  not_defor_counterfactual <- rowSums(year_counterfactual)<1 *1
  defor_year_counterfactual <- max.col(year_counterfactual, ties.method = "first") 
  defor_counterfactual <- transform(year_counterfactual, defor_year_counterfactual = ifelse(not_defor_counterfactual==1, years*2+1, defor_year_counterfactual))
  defor_counterfactual <- tibble::rownames_to_column(defor_counterfactual)
  names(defor_counterfactual)[1] <- paste("pixels")
  
  panels_counterfactual <- defor_counterfactual %>%
    dplyr::select(pixels, defor_year_counterfactual) %>%
    inner_join(panels_counterfactual, by = "pixels")
  
  cols.num <- c("pixels", "grid", "property", "county", "year")
  panels_counterfactual[cols.num] <- sapply(panels_counterfactual[cols.num],as.numeric)
  
  panels_counterfactual <- panels_counterfactual %>%
    mutate(indic = year - defor_year_counterfactual) %>%
    mutate(defor_counter = ifelse(indic > 0, 1, y_counterfactual)) %>%
    mutate(y_it = ifelse(indic > 0, NA, y_counterfactual) )  
  
  #panels_counterfactual$defor_counter <- as.factor(panels_counterfactual$defor_counter)
  
  panel_extra <- panels %>%
    dplyr::select(pixels, year, defor)
  
  extra_defor <- panels_counterfactual %>%
    dplyr::select(pixels, year, treat, defor_counter, geometry) %>%
    inner_join(panel_extra, by = c("pixels", "year")) 
  
  extra_defor$defor_counter <- as.numeric(as.character(extra_defor$defor_counter))
  extra_defor$defor<- as.numeric(as.character(extra_defor$defor))
  
  ######################################################################
  ################# Creating Landscape Maps  ###########################
  ######################################################################
  
  
  fills <- c("pixels not deforested in intervention area" = "#009E73", 
             "pixels not deforested in untreated area" = "#56B4E9", 
             "deforested pixels" = "#F0E442", 
             "deforestation in first period" = "red", 
             "deforested pixels in couterfactual without intervention" = "#D55E00",
             "deforestation in second period" = "gray30",
             "pixels not deforested" = "#009E73")
  
  shapes <- c("deforested pixel" = 22)
  
  colors <- c("county boundaries" = "black", 
              "property boundaries" = "darkgoldenrod4"
              #"deforested pixels in couterfactual without intervention" = "red"
              )
  
  panels$defor <- as.factor(panels$defor)
  
  # plot_df_period1 <- panels %>%
  #   st_as_sf() %>%
  #   dplyr::select(pixels, year, treat, defor) %>%
  #   filter(year == 1  & defor == 1)
  # 
  # landscape_period1 <- 
  #   ggplot() + 
  #   geom_sf(data = intervention_area, aes(fill = "intervention area"), color = "#a1d76a") +
  #   geom_sf(data = control_area, aes(fill = "control area"), color = "lightblue")+
  #   geom_sf(data = plot_df_period1, aes(fill = "deforested pixel"), color = "NA", shape = 22, alpha = .9, size = 1.5)+
  #   geom_sf(data = p_bounds, aes(color = "property boundaries"), fill = "NA")+
  #   geom_sf(data = c_bounds, aes(color = "county boundaries"), size = 1, fill = "NA")+
  #   scale_color_manual(values = colors) +
  #   scale_fill_manual(values= fills) +
  #   theme(panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank(), 
  #         panel.background = element_rect(fill = "gray95"), 
  #         legend.title = element_blank(),
  #         axis.text = element_blank(),
  #         axis.ticks = element_blank(),
  #         legend.position = "none")
  # 
  plot_df_period2 <- panels %>%
    st_as_sf() %>%
    dplyr::select(pixels, year, treat, defor, y_it) %>%
    filter(year == max(year)  & defor == 1)
  
  landscape_period2 <- 
    ggplot() + 
    geom_sf(data = intervention_area, aes(fill = "pixels not deforested in intervention area"), color = "#a1d76a") +
    geom_sf(data = control_area, aes(fill = "pixels not deforested in untreated area"), color = "lightblue")+
    #geom_sf(data = plot_df_period1, aes(fill = "deforestation in first period"), color = "NA", shape = 22, alpha = .9, size = 1)+
    geom_sf(data = plot_df_period2, aes(fill = "deforested pixels"), color = "NA", shape = 22, alpha = .9, size = 1.2)+
    geom_sf(data = p_bounds, aes(color = "property boundaries"), fill = "NA")+
    geom_sf(data = c_bounds, aes(color = "county boundaries"), size = 1, fill = "NA")+
    scale_fill_manual(values= fills) +
    scale_color_manual(values = colors) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white"), 
          legend.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key=element_blank())
  
  # landscape_period2_nolegend <- landscape_period2 +
  #   theme(panel.grid.major = element_blank(), 
  #         panel.grid.minor = element_blank(), 
  #         panel.background = element_rect(fill = "gray95"), 
  #         legend.title = element_blank(),
  #         axis.text = element_blank(),
  #         axis.ticks = element_blank(),
  #         legend.position = "none")
  
  extra_period2 <- extra_defor %>%
    mutate(extra = defor_counter - defor) %>%
    st_as_sf() %>%
    filter(year == max(year)  & extra == 1 &treat==1)
  extra_period2$extra<- as.factor(extra_period2$extra)
  
  landscape_period2_counter <- landscape_period2 +
    geom_sf(data = extra_period2, aes(fill = "deforested pixels in couterfactual without intervention"), color = "NA", shape = 22, alpha = .9, size = 1.2)+
    guides(shape=FALSE)
  
  ######################################################################
  outputs = list("landscape_period1" = landscape_period1, "landscape_period2" = landscape_period2, "landscape_period2_counter" = landscape_period2_counter)
  return(outputs)
  
  #end function  
}  
