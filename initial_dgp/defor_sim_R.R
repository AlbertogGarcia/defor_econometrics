# transitioning initial python simulation to R
#created 14 may 2019

#Simulation model to test for possible bias in generalized DID model (two way FE) w/ binary outcome, dropped observations

#packages
library(dplyr)
library(purrr)
library(tidyverse)
library(reshape2)
library(data.table)
library(sandwich)
library(plm)
library(lmtest)
library(clubSandwich)

#years <- 2               # number of years in each of two periods
#nobs <- 10000           # number of observations in each of two groups

#starting function
defor_sim <- function(nobs, years){


#define parameters
          
d_outside <- 0.30        # true deforestation rate outside treated area in period 2 
diff <- 0.40             # pre-treatment difference between treatment and control
trend <- -0.10           # trend in deforestation rate across periods
att <- -0.16             # average treatement effect on the treated


#deforestation rates
m_00 <- d_outside                        # average deforestation rate outside treated area in p1
m_01 <- d_outside + trend                # average deforestation rate outside treated area in p2
m_10 <- d_outside + diff                 # average deforestation rate inside treated area in p1
m_11 <- d_outside + diff + trend + att   # average deforestation rate inside treated area in p2


# introducing three types of randomness
#     i: Individual-level
#     y: year-level changes
#     e: common structure across all observations

std_i <- 0.1
std_y <- 0.01
std_e <- 0.01

##########################  This section could probably be more efficiently coded  ###############################
i_err <- rnorm(nobs*2 , 0 , std_i)
y_err <- rnorm(years*2 , 0 , std_y)
e_err <- matrix(rnorm(nobs*2 * years*2 , 0 , std_e), ncol = years*2)

i_err <- replicate(years*2, i_err)
y_err <- t(replicate(nobs*2, y_err))

df <- data.frame( i_err + y_err + e_err)
###################################################################################################################

#Add in average deforestation rates for each group
df[1:nobs,1:years] <- df[1:nobs,1:years] + m_00
df[1:nobs,(years+1):ncol(df)] <- df[1:nobs,(years+1):ncol(df)] + m_01
df[(nobs+1):nrow(df),1:years] <- df[(nobs+1):nrow(df),1:years] + m_10
df[(nobs+1):nrow(df),(years+1):ncol(df)] <- df[(nobs+1):nrow(df),(years+1):ncol(df)] + m_11


#define treatment
df$treat <- rep(0, nrow(df))
df$treat[(nobs+1):nrow(df)] <- df$treat[(nobs+1):nrow(df)]+1

df <- tibble::rownames_to_column(df)
df <- unite(df, idx, treat, rowname, sep = "_", remove = FALSE)

df <- unite(df, idx, idx, treat, sep = " , ", remove = TRUE)

rownames(df) <- df$idx
df <- subset(df, select = -c(rowname,idx))

### simulating deforestation ###
defor_draw <- matrix(runif(nobs*2 * years*2 , 0 , 1), ncol = years*2)
defor_df <- data.frame( df > defor_draw)*1
not_defor <- rowSums(defor_df)<1 *1
defor_year <- max.col(defor_df, ties.method = "first")       #creating defor_year variable


defor_df <- cbind(defor_df, defor_year)
defor_df <- transform(defor_df, defor_year = ifelse(not_defor==1, years*2+1, defor_year))

#names(defor_df)[1:(years*2)] <- paste("defor", colnames(defor_df[1:(years*2)]), sep = "")

defor_df <- tibble::rownames_to_column(defor_df)

defor_df <- separate(defor_df, rowname, c("idx", "treat"), sep= ",", remove=TRUE)

names(defor_df)[3:(years*2+2)] <- c(1:(years*2))

defor_df <- melt(defor_df, id.vars = c('idx', 'treat', 'defor_year'))


setnames(defor_df, old=c("variable","value"), new=c("year", "defor"))

#replace defor with NA for any year > defor_year
defor_df$year <- as.numeric(defor_df$year)
defor_df$indic <- (defor_df$year - defor_df$defor_year)
defor_df$defor <- ifelse(defor_df$indic > 0 , NA, defor_df$defor)
defor_df <- subset(defor_df, select = -c(indic))


defor_df$post <- (defor_df$year > years)*1

assign('defor_df',defor_df, envir=.GlobalEnv)
}





#defor_df <- unite(defor_df, index, year, idx, sep = ",", remove = TRUE)
#rownames(defor_df) <- defor_df$index
#defor_df <- subset(defor_df, select = -c(index))

# DID <- lm(defor ~  post*treat, 
#           data = defor_df
# )
# coeftest(DID, vcov. = vcovCL)
# 
# twoway.fe <- plm(defor ~  post*treat, 
#              data   = defor_df, 
#              method = "within", #fixed effects model
#              effect = "twoway", #unit and year fixed effects
#              index  = c("idx", "year")
# )

