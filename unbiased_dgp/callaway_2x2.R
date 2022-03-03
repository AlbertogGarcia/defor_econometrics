# this function's intention is to provide distributional parameters for the various specifications that have the typical binary outcome variable
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(fixest)
library(did)
source(here::here('unbiased_dgp', 'deforestation_DGP.R'))
#source('deforestation_DGP.R')
#begin function

callaway_2x2 <- function(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a=.1, std_v=.25, std_p=0, cellsize, ppoints, cpoints){
  
  countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
  pixloc_df = countyscape$pixloc_df
  ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)
  
  #preallocate n x 4 matrix
  coeffmatrix <- matrix(nrow = n, ncol = 4)
  
  n_mod = 4
  summ_row <- n_mod * n
  
  summary_long <- data.frame('b0'= rep(b0, summ_row), 'b1'= rep(b1, summ_row), 'b2_0'= rep(b2_0, summ_row), 'b2_1'= rep(b2_1, summ_row), 'b3'= rep(b3, summ_row), 
                             'std_a'= rep(std_a, summ_row), 'std_v'= rep(std_v, summ_row), 'std_p'= rep(std_p, summ_row),
                             'iteration' = rep(NA, summ_row), 
                             'model'=rep(NA, summ_row),
                             'bias'=rep(NA, summ_row), 
                             stringsAsFactors=FALSE)
  
  for(i in 1:n){
    Nobs <- length(pixloc$treat)  
    panels <- fabricate(
      pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
      year = add_level(N = (years*2), nest = FALSE),
      obs = cross_levels(
        by = join(pixels, year),
        post = ifelse(year > years, 1, 0),
        treat.year = ifelse(treat==1, years+1, 0),
        v_it = rnorm(N, 0, std_v),
        ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it
      )
    )
    
    #generate random 
    errortable_property <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, std_p))
    #errortable_county <- data.frame(county = as.character(unique(pixloc$county)), c_err = rnorm(length(unique(pixloc$county)), 0, std_c))
    
    panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)
    
    panels <- panels %>%
      inner_join(pixloc, by = c("pixels", "treat")) %>%
      inner_join(errortable_property, by = "property") %>%
      #inner_join(errortable_county, by = "county") %>%
      mutate(ystar = ystar + p_err) %>%
      mutate(y = (ystar > 0)*1 )
    
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
      mutate(indic = year - defor_year) %>%
      mutate(defor = ifelse(indic > 0, 1, y))%>%
      mutate(y_it = ifelse(indic > 0, NA, y))
    
    # aggregate up to county in each year 
    suppressWarnings(
      countylevel_df <-  aggregate(panels, by = list(panels$county, panels$treat, panels$year, panels$treat.year), FUN = mean, drop = TRUE)[c("county", "property", "treat", "post", "year","defor", "carea", "treat.year")]
    )
    
    
    countylevel_df <- countylevel_df[order(countylevel_df$county, countylevel_df$year),]
    countylevel_df <- slide(countylevel_df, Var = "defor", GroupVar = "county", NewVar = "deforlag",
                            slideBy = -1, reminder = FALSE)
    
    ##### creating forested share variable #####
    countylevel_df$forshare <- 1 -  countylevel_df$defor
    countylevel_df$forsharelag <- 1 -  countylevel_df$deforlag
    
    #generate outcome var
    countylevel_df$deforrate <- ((countylevel_df$forsharelag- countylevel_df$forshare) / countylevel_df$forsharelag)
    #remove any infinite values
    
    countylevel_df <- 
      countylevel_df %>% 
      filter_all(all_vars(!is.infinite(.)))%>%
      drop_na(deforrate)
    # call defor_sim function to simulate dataframe, returned as panels  
    
    
    coeffmatrix[i,1]  <- feols(y_it ~  post:treat|year+pixels, data = panels
    )$coefficients - ATT
    
    # run two-way fixed effects    
    coeffmatrix[i, 2] <- feols(deforrate ~  post:treat|year+county, data = countylevel_df
    )$coefficients - ATT
    
    pixel_attgt <- att_gt(yname = "y_it",
                           tname = "year",
                           idname = "pixels",
                           gname = "treat.year",
                           control_group = "nevertreated",
                           data = panels,
                           bstrap=FALSE)
    
    county_attgt <- att_gt(yname = "deforrate",
                          tname = "year",
                          idname = "county",
                          gname = "treat.year",
                          control_group = "nevertreated",
                          data = countylevel_df,
                          bstrap=FALSE)
    
    coeffmatrix[i, 3] <- aggte(pixel_attgt)$overall.att - ATT
    
    coeffmatrix[i, 4] <- aggte(county_attgt)$overall.att - ATT
      
    
    firstcol = which(colnames(summary_long)=="iteration")
    lastcol = which(colnames(summary_long)=="bias")
    
    summary_long[i,c(firstcol:lastcol)] <- c(
      i,
      "TWFE_pixel",
      coeffmatrix[i,1]
    )

    summary_long[i+n,c(firstcol:lastcol)] <- c(
      i,
      "TWFE_county",
      coeffmatrix[i,2]
    )

    summary_long[i+n*2,c(firstcol:lastcol)] <- c(
      i,
      "callaway_pixel",
      coeffmatrix[i,3]
    )

    summary_long[i+n*3,c(firstcol:lastcol)] <- c(
      i,
      "callaway_county",
      coeffmatrix[i,4]
    )
    
    
    #end for loop
    print(i)
  }  
  
  # get distribution information from matrix  
  
  b_coeff <- as.data.frame(coeffmatrix)
  actual <- rep(0, times = n)
  
  names(b_coeff)[1] <- paste("TWFE_pixel")
  names(b_coeff)[2] <- paste("TWFE_county")
  names(b_coeff)[3] <- paste("callaway_pixel")
  names(b_coeff)[4] <- paste("callaway_county")
  
  suppressWarnings(cbias <- melt(b_coeff, value.name = "bias"))
  
  summary_long <- summary_long %>%
    select(iteration, everything())
  
  
  outputs = list("did_biases" = b_coeff, "summary_long" = summary_long)
  return(outputs)
  
  #end function  
}  




