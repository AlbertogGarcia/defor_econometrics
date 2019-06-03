# econometric analyses
# 

library(sandwich)
library(plm)
library(lmtest)
library(clubSandwich)
library(DataCombine)


################ dropping observations in periods t+1...T when defor_it = 1 #########
#defor_sim(nobs, years)
defor_sim(10000, 3) #initial simulation with 10000 observations in each group/2 periods in pre/post
# spits out dataframe defor_df


#simple DID 
DID_drop <- lm(defor ~  post*treat, 
          data = defor_df
)


#two-way fixed effects model
twoway.fe <- plm(defor ~  post*treat, 
                 data   = defor_df, 
                 method = "within", #fixed effects model
                 effect = "twoway", #unit and year fixed effects
                 index  = c("idx", "year")
)

################ setting observations = 1 in periods t+1...T when defor_it = 1 #########

defor_df$defor[is.na(defor_df$defor)] <- 1







                


################ with property level perturbations ########################

prop_deforsim(10000, 3, 500, 20) #initial simulation with 10000 observations in each group/2 periods in pre/post
# spits out both county level and property level dataframes countypert_df and proppert_df

year1 <- subset(countyprop_df, year== 1)
colnames(year1)[colnames(year1)=="defor"] <- "defor0"
year1 <- subset(year1, select = c(countyid, defor0))
st_geometry(year1) <- NULL
countyprop_df <-  merge(countyprop_df, year1, by = "countyid")
countyprop_df$deforrate1 <- countyprop_df$defor / countyprop_df$defor0

year1 <- subset(proppert_df, year== 1)
colnames(year1)[colnames(year1)=="defor"] <- "defor0"
year1 <- subset(year1, select = c(property, defor0))
st_geometry(year1) <- NULL
proppert_df <-  merge(proppert_df, year1, by = "property")

# proppert_df$defor0 <- ifelse(proppert_df$defor0 >0, proppert_df$defor0, .1 )
# 
# proppert_df$deforrate1 <- proppert_df$defor / proppert_df$defor0


DID_county <- lm(deforrate1 ~  post*treat, 
               data = countyprop_df
)

summary(DID_county)


# DID_prop <- lm(deforrate1 ~  post*treat, 
#                data = proppert_df
# )
# 
# 
# summary(DID_prop)




#two-way fixed effects model
twoway.fec <- plm(deforrate1 ~  post*treat, 
                 data   = countyprop_df, 
                 method = "within", #fixed effects model
                 effect = "twoway", #unit and year fixed effects
                 index  = c("countyid", "year")
)

summary(twoway.fec)
# 
# twoway.fep <- plm(deforrate1 ~  post*treat, 
#                   data   = proppert_df, 
#                   method = "within", #fixed effects model
#                   effect = "twoway", #unit and year fixed effects
#                   index  = c("property", "year")
# )
# 
# summary(twoway.fep)
