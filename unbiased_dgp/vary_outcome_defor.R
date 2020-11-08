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
source('grid_landscape.R')

vary_outcome_defor <- function(n, nobs, years, ATT, trend, base_start, fixed_diff = 0,increments = 0.05, steps = 10, std_a = 0.1, std_v = 0.25, cellsize){
  
  base_range <- seq(from = base_start, to = (base_start + (increments*steps)), by = increments)
  base_list <- as.list(base_range)
  
  #preallocate n x ___ matrix
  out_matrix1 <- matrix(nrow = n, ncol = length(base_list))
  out_matrix2 <- matrix(nrow = n, ncol = length(base_list))
  out_matrix3 <- matrix(nrow = n, ncol = length(base_list))
  
  gridscape = grid_landscape(nobs, cellsize)
  pixloc_df = gridscape$pixloc_df
  gridcoords = gridscape$gridcoords
  pixloc <- pixloc_df
  
  for(k in base_list){
    tic("loop")
    base_0 <- k
    base_1 <- k + fixed_diff
    place <- which(base_list == k)
    
    #compute relevant beta parameters
    std_av = (std_a^2+std_v^2)^.5
    b0 = qnorm(base_0, mean = 0, sd = std_av)
    b1 = qnorm(base_1, mean = 0, sd = std_av) - b0
    b2_0 = qnorm(trend + base_0, mean = 0, sd = std_av) - b0
    b2_1 = qnorm(trend + base_1, mean = 0, sd = std_av) - b0 - b1
    b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_av) + ATT , mean = 0, sd = std_av) - (b0 + b1 + b2_1)
    
    for(i in 1:n){
      
      Nobs <- length(pixloc$treat)   
      panels <- fabricate(
        pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
        year = add_level(N = (years*2), nest = FALSE),
        obs = cross_levels(
          by = join(pixels, year),
          post = ifelse(year > years, 1, 0),
          v_it = rnorm(N, 0, std_v),
          ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
          y = (ystar > 0)*1
        )
      )
      
      panels$year <- as.numeric(panels$year)
      panels <- panels %>%
        mutate(post = ifelse(year > years, 1, 0))
      
      #need to determine which year deforestation occurred
      year_df <- subset(panels, select = c(pixels, year, y))
      #year_df <- melt(year_df, id.vars = c("pixels", "y_it"), value.name = "year")
      year_df <- dcast(year_df, pixels ~ year , value.var = "y")
      rownames(year_df) <- year_df$pixels
      year_df <- subset(year_df, select = -c(pixels))
      
      #creating variable for the year a pixel is deforested
      not_defor <- rowSums(year_df)<1 *1
      defor_year <- max.col(year_df, ties.method = "first") 
      defor_df <- transform(year_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))
      defor_df <- tibble::rownames_to_column(defor_df)
      names(defor_df)[1] <- paste("pixels")
      defor_df <- subset(defor_df, select = c(pixels, defor_year))
      panels <- merge(defor_df, panels, by = "pixels")
      
      
      # creating three outcome variables for each possible situation
      ### y: allows the outcome to switch between 0 and 1 across years
      ### y_it: outcome is dropped in years after pixel is first deforested
      ### defor: outcome is set to 1 in each year after the pixel is deforested
      panels$year <- as.numeric(panels$year)
      panels$indic <- (panels$year - panels$defor_year)
      panels$defor <- ifelse(panels$indic > 0 , 1, panels$y)
      panels <- subset(panels, select = -c(indic))
      
      panels$pixels <- as.numeric(gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE))
      
      
      panels <- panels %>%
        inner_join(pixloc, by = c("pixels", "treat"))
      
      
      # aggregate up to county in each year 
      suppressWarnings(
        gridlevel_df <-  aggregate(panels, by = list(panels$grid, panels$treat, panels$year), FUN = mean, drop = TRUE)[c("grid", "treat", "post", "year","defor")]
      )
      
      gridlevel_df <- gridlevel_df[order(gridlevel_df$grid, gridlevel_df$year),]
      
      
      gridlevel_df <- gridlevel_df %>%
        slide(Var = "defor", GroupVar = "grid", NewVar = "deforlag", slideBy = -1, reminder = FALSE) %>%
        mutate(forshare = (1-defor)) %>%
        mutate(forsharelag = (1-deforlag))
      
      #### create baseline defor and forshare vars 
      #county level
      year1 <- subset(gridlevel_df, year== 1)
      colnames(year1)[colnames(year1)=="defor"] <- "defor0"
      year1 <- subset(year1, select = c(grid, defor0))
      #st_geometry(year1) <- NULL
      gridlevel_df <-  gridlevel_df %>%
        inner_join( year1, by = "grid") %>%
        mutate(forshare0 = (1 - defor0))
      
      #generate outcome var
      
      gridlevel_df <- gridlevel_df %>%
        mutate(deforrate1 = (forsharelag - forshare) / forsharelag) %>%
        mutate(deforrate2 = (forsharelag - forshare) / forshare0) %>%
        mutate(deforrate3 = log(forsharelag / forshare))
      
      
      #remove any infinite values
      #gridlevel_df <- subset(gridlevel_df, select = -c(geometry))
      gridlevel_df <- gridlevel_df %>% 
        filter_all(all_vars(!is.infinite(.)))
      
      # run two-way fixed effects with outcome 1 
      out_matrix1[i,place] <- plm(deforrate1 ~  post*treat, 
                              data   = gridlevel_df, 
                              method = "within", #fixed effects model
                              effect = "twoway", #grid and year fixed effects
                              index  = c("grid", "year")
      )$coefficients - ATT
      
      out_matrix2[i,place] <- plm(deforrate2 ~  post*treat, 
                              data   = gridlevel_df, 
                              method = "within", #fixed effects model
                              effect = "twoway", #grid and year fixed effects
                              index  = c("grid", "year")
      )$coefficients - ATT
      
      out_matrix3[i,place] <- plm(deforrate3 ~  post*treat, 
                              data   = gridlevel_df, 
                              method = "within", #fixed effects model
                              effect = "twoway", #grid and year fixed effects
                              index  = c("grid", "year")
      )$coefficients - ATT
      
    } # end of n loop
    
    print(k)
    toc()
    
  } # end of list loop
  
  outbias_df1 <- as.data.frame(cbind(colMeans(out_matrix1), base_range)) 
  outbias_df2 <- as.data.frame(cbind(colMeans(out_matrix2), base_range)) 
  outbias_df3 <- as.data.frame(cbind(colMeans(out_matrix3), base_range)) 
  
  colors <- c("outcome 1" = "red", 
              "outcome 2" = "blue",
              "outcome 3" = "green")
  
  plot <- ggplot() + 
    geom_line(data = outbias_df1, aes(x = base_range, y = V1, color="outcome 1"), size =1) +
    geom_line(data = outbias_df2, aes(x = base_range, y = V1, color="outcome 2"), size =1) +
    geom_line(data = outbias_df3, aes(x = base_range, y = V1, color="outcome 3"), size =1) +
    labs(x = "pre-treatment deforestation rates", y = "Bias", caption = paste("Bias depends on deforestation formula")) +
    scale_color_manual(values = colors)+
    geom_hline(yintercept = 0, linetype = "dashed")
  
  outputs = list("outbias_df1" = outbias_df1, "outbias_df2" = outbias_df2, "outbias_df3" = outbias_df3, "plot" = plot)
  
  
} # end of function
    