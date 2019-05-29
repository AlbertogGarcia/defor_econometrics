# econometric analyses
# may want to incorporate monte carlos to see distribution of point estimates at some point

library(sandwich)
library(plm)
library(lmtest)
library(clubSandwich)

################ dropping observations in periods t+1...T when defor_it = 1 #########

defor_sim(nobs, years) #initial simulation with nobs observations in each group/years periods in pre/post
# spits out dataframe defor_df


#simple DID 
DID_pixel <- lm(defor ~  post*treat, 
          data = defor_df
)

#two-way fixed effects model
twoway.fe <- plm(defor ~  post*treat, 
                 data   = defor_df, 
                 method = "within", #fixed effects model
                 effect = "twoway", #unit and year fixed effects
                 index  = c("idx", "year")
)



################ aggregation without property or county level perturbations ########################
####### defor set to 1 after first occurence of deforestation

grid_fcn(defor_df, cellsize) 
#spits out grid level (counties)dataframe countylevel_df

property_fcn(defor_df, cellsize, psize)
#spits out property level dataframe without perturbations propertylevel_df





                


################ with property level perturbations ########################

prop_deforsim(nobs, years, cellsize, psize) #initial simulation with nobs observations in each group/years periods in pre/post
# spits out both county level and property level dataframes countypert_df and proppert_df




