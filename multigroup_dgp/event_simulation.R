#begin function
library(tidyverse)
library(tictoc)
library(fixest)
library(here)
library(DeclareDesign)
library(did2s)

event_simulation <- function(n, pixloc, st_a, std_v, std_p, b0a, b1a, b2a, b3a, b4a, b0b, b1b, b2b, b3b, b4b, b0c, b1c, b2c, b3c, b4c, tau_a, tau_a2, tau_a3, tau_b, tau_b2){
  
  pixel_es_long <- data.frame(matrix(ncol = 5, nrow = 0))
  county_es_long <- data.frame(matrix(ncol = 5, nrow = 0))
  
  #provide column names
  cnames <- c('term', 'estimate', 'std.error', 'estimator', 'iteration')
  colnames(pixel_es_long) <- cnames
  colnames(county_es_long) <- cnames
  
  
  for(i in 1:n){
    
    Nobs <- length(pixloc$G)  
    
    panels <- fabricate(
      pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), G = pixloc$G),
      year = add_level(N = 5, nest = FALSE),
      obs = cross_levels(
        by = join(pixels, year),
        post = (year >= G)*1,
        GU = (G==0)*1,
        G3 = (G==3)*1,
        G4 = (G==4)*1,
        v_it = rnorm(N, 0, std_v),
        b0 = b0a*G3 + b0b*G4 + b0c*GU,
        b1 = b1a*G3 + b1b*G4 + b1c*GU,
        b2 = b2a*G3 + b2b*G4 + b2c*GU,
        b3 = b3a*G3 + b3b*G4 + b3c*GU,
        b4 = b4a*G3 + b4b*G4 + b4c*GU,
        tau = tau_a*G3*post + tau_b*G4*post ,
        tau_2 = tau_a2*G3*((year >=G+1)*1) + tau_b2*G4*((year >=G+1)*1),
        tau_3 = tau_a3*G3*((year >=G+2)*1),
        ystar = b0 + b1*((year >=2)*1) + b2*((year >=3)*1) + b3*((year >=4)*1) + b4*((year >=5)*1) + tau +tau_2 + tau_3 + a_i + v_it,
        y = (ystar > 0)*1
      )
    )
    
    #generate random 
    errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "G")) %>%
      inner_join(errortable_property, by = "property") %>%
      #inner_join(errortable_county, by = "county") %>%
      mutate(ystar = ystar + p_err) %>%
      mutate(y = (ystar > 0)*1,
             pixels = as.numeric(pixels),
             year = as.numeric(year),
             defor_indic = ifelse(y==1, year, 99))%>%
      group_by(pixels)%>%
      mutate(defor_year = min(defor_indic),
             defor = (year>=defor_year)*1,
             y_it = ifelse(year>defor_year, NA, defor)
      )
    
    
    
    #########################################################################
    ######### pixel estimates  
    #########################################################################
    
    y_it_es <- my_event_study(yname = "y_it",
                              tname = "year",
                              idname = "pixels",
                              gname = "G",
                              data = panels) 
    
    pixel_es_long <- rbind(pixel_es_long, y_it_es)%>%
      mutate(std.error = as.numeric(ifelse(is.na(std.error), 0.0, std.error)),
             estimate = as.numeric(estimate),
             term = as.numeric(term),
             uoa = "pixel")
    
    #plot_event_study(pixel_es_long, seperate = TRUE, horizon = NULL)
    
    #########################################################################
    ######### county level dataframe  
    #########################################################################
    
    
    countylevel_df <- panels %>%
      group_by(county, year, G) %>% 
      dplyr::summarise(defor = mean(defor),
                       G = mean(G),
                       county = as.numeric(county))%>%
      distinct()
    
    countylevel_df <- countylevel_df[order(countylevel_df$county, countylevel_df$year),]
    countylevel_df <- slide(countylevel_df, Var = "defor", GroupVar = "county", NewVar = "deforlag",
                            slideBy = -1, reminder = FALSE)%>%
      mutate(forshare = 1-defor,
             forsharelag = 1 - deforlag,
             deforrate = (forsharelag- forshare) / forsharelag
      ) %>% 
      filter_all(all_vars(!is.infinite(.)))%>%
      drop_na(deforrate)
    
    #########################################################################
    ######### county estimates  
    #########################################################################
    
    county_es <- my_event_study(yname = "deforrate",
                                tname = "year",
                                idname = "county",
                                gname = "G",
                                data = countylevel_df) 
    
    county_es_long <- rbind(county_es_long, county_es)%>%
      mutate(std.error = as.numeric(ifelse(is.na(std.error), 0.0, std.error)),
             estimate = as.numeric(estimate),
             term = as.numeric(term),
             uoa = "county")
    
    #plot_event_study(county_es_long, seperate = TRUE, horizon = NULL)
    
  }  
  
  es_long <- rbind(county_es_long, pixel_es_long)
  
  es_long[] <- lapply(es_long, as.character)
  
  return("es_long" = es_long)
  
}
