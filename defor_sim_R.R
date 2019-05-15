# transitioning initial python simulation to R
#created 14 may 2019

#Simulation model to test for possible bias in generalized DID model (two way FE) w/ binary outcome, dropped observations

#packages
library(dplyr)
library(purrr)
library(tidyverse)


#define parameters
years <- 2               # number of years in each of two periods
nobs <- 10000            # number of observations in each of two groups
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
e_err <- matrix(rnorm(nobs*2 * years*2 , 0 , std_e), ncol = 4)

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


defor_draw <- matrix(runif(nobs*2 * years*2 , 0 , 1), ncol = 4)
defor_df <- data.frame( df > defor_draw)*1



