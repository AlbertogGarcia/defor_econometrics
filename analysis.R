# econometric analyses
# 

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
####### defor set to 1 after first occurence of deforestation (nothing dropped)


property_fcn(defor_df, cellsize, psize)
#spits out property and county level dataframe without perturbations propertylevel_df and countylevel_df

deforrate1 <-                            #defor_t/defor_{t-1}
deforrate2 <-                            #defor_t/defor_0  
  



                


################ with property level perturbations ########################

prop_deforsim(nobs, years, cellsize, psize) #initial simulation with nobs observations in each group/years periods in pre/post
# spits out both county level and property level dataframes countypert_df and proppert_df




