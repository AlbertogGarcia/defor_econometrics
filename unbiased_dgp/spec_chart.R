#####################################################################
#### r script that generates specification charts for paper
#####################################################################
library(dplyr)
library(Metrics)
source(here::here('unbiased_dgp', 'schart.R'))
summary_pweights <- readRDS("unbiased_dgp/summary_pweights.rds")
specchart_long <- summary_pweights

summary <- specchart_long %>%
  mutate_at(vars(bias, cover, ATT_diff, ls_ATT, p_ATT), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes) & pixel.fe == FALSE)%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover),
            ls_ATT = mean(ls_ATT),
            ATT_diff = mean(ATT_diff),
            p_ATT = mean(p_ATT))%>%
  select(Bias, everything())

df_summary <- summary %>%
  select(-c(ls_ATT, p_ATT, ATT_diff))

#df_summary <- rbind(df_summary[11,], df_summary[4:10,], df_summary[1:3,])

par(oma=c(1,0,1,1))

labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "Weights:" = c("area weights"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))
f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))
coverage <- f(round(df_summary$cover, digits=3))
RMSE <- f(round(df_summary$RMSE, digits=5))
test <- as.data.frame(subset(df_summary, select=-c(cover, RMSE)))

test <- subset(test, select=-c(se_pixel, se_grid, se_property, se_county))%>%
  distinct()

labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "Weights:" = c("area weights")
)

index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart
ylim <- c(-0.0125,0.01)
#bottomline = (min(ylim)+topline)/2
#Create the plot

schart(test,labels, ylim = ylim, index.ci=index.ci#,# col.est=c("grey50","#00A1D5")
       , ylab="Bias                              "
       , highlight=c(3,4), 
       col.est=c("black","#00A1D5"),
       col.dot=c("black","lightgrey","red","#00A1D5"),
       bg.dot=c("black","lightgrey","white","#00A1D5")
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
abline(h=mean(summary$ATT_diff), col = "#00A1D5", lty = "dashed")
mtext("property level\nATT", side=2, at = mean(summary$ATT_diff), font=2, las=1, line=.5, cex = .75)
mtext("landscape\nATT", side=2, at = 0, font=2, las=1, line=.5, cex = .75)



#####################################################################
#### initial spec chart
#####################################################################


specchart_long <- readRDS("unbiased_dgp/summary_long.rds")

df_summary <- specchart_long %>%
  mutate_at(vars(bias, cover), as.numeric
         )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes))%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())

df_summary <- rbind(df_summary[11,], df_summary[4:10,], df_summary[1:3,])

par(oma=c(1,0,1,1))

labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))
f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))
coverage <- f(round(df_summary$cover, digits=3))
RMSE <- f(round(df_summary$RMSE, digits=5))
test <- as.data.frame(subset(df_summary, select=-c(cover, RMSE)))


index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart
ylim <- c(-0.035,0.04)

topline = -0.015

midline = topline-0.02-.0025

ylim <- c(midline-.02,0.04)
#bottomline = (min(ylim)+topline)/2
#Create the plot

schart(test,labels, ylim = ylim, index.ci=index.ci, col.est=c("grey50","#D55E00")
       , ylab="Bias", highlight=1
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c( 0, 0.03), labels=TRUE)
abline(h=topline)
abline(h=midline)
#abline(h=bottomline)
lapply(1:length(RMSE), function(i) {
  # text(x= i, y=min(ylim)+.002, paste0(i), col="black", font=2, cex = 1)
  mtext(paste0(i), side=1, at = i, font=2, cex=.95)#, line=1, at=-1)
  text(x= i, y=midline+0.005, paste0(RMSE[i]), col="black", font=1, cex=.9)
  text(x= i, y=min(ylim)+0.0025, paste0(coverage[i]), col="black", font=1, cex=.9 )
})
# mtext('coverage\nprobability', side=2, las=1, at = loc_1, adj = 1, font=2, line=1)
# mtext('RMSE', side=2, las=1, at = loc_2, font=2, line=1)
text(x=mean(1:nrow(test))
     , y=midline-.005, "coverage probability (95% desired)", col="black", font=2)
text(x=mean(1:nrow(test))
     , y=topline-.005, "RMSE", col="black", font=2)
legend(x=6.5, y=0.05, col = "grey50", legend = "0.05 to 0.95 quantile \n of estimate distribution", seg.len=0.65, inset = 0.005,  box.lty=0, cex=0.8, lty = 1, lwd = 4, bg="transparent")



###############################################################################################
###### spec chart to highlight selection bias
###############################################################################################

specchart_long <- readRDS("unbiased_dgp/summary_selection.rds")

df_summary <- specchart_long %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe,grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes))%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe,property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())%>%
  filter(pixel.fe ==FALSE)

#df_summary <- rbind(df_summary[11,], df_summary[4:10,], df_summary[1:3,])

par(oma=c(1,0,1,1))

labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c( "grid FE", "property FE", "county FE", "treatment FE"),
               "Weights:" = c("area"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))
f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))
coverage <- f(round(df_summary$cover, digits=3))
RMSE <- f(round(df_summary$RMSE, digits=5))
test <- as.data.frame(subset(df_summary, select=-c( pixel.fe, cover, RMSE)))

highlight_rows <- which(test[ , "pixel"] == TRUE)

index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart
ylim <- c(-0.035,0.015)

topline = -0.015

midline = topline-0.01-.0025

ylim <- c(midline-.01,0.015)
#bottomline = (min(ylim)+topline)/2
#Create the plot

schart(test,labels, ylim = ylim, index.ci=index.ci, 
       col.est=c("grey50","#DF8F44"),
       col.dot=c("grey50","lightgrey","red","#DF8F44"),
       bg.dot=c("grey50","lightgrey","white","#DF8F44")
       , ylab="Mean bias and\n.05 to .95 quantiles", highlight=highlight_rows
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE)
abline(h=topline)
abline(h=midline)
#abline(h=mean(df_summary$sel_bias))
#abline(h=bottomline)
lapply(1:length(RMSE), function(i) {
  # text(x= i, y=min(ylim)+.002, paste0(i), col="black", font=2, cex = 1)
  mtext(paste0(i), side=1, at = i, font=2, cex=.95)#, line=1, at=-1)
  text(x= i, y=midline+0.005, paste0(RMSE[i]), col="black", font=1, cex=.9)
  text(x= i, y=min(ylim)+0.0025, paste0(coverage[i]), col="black", font=1, cex=.9 )
})
# mtext('coverage\nprobability', side=2, las=1, at = loc_1, adj = 1, font=2, line=1)
# mtext('RMSE', side=2, las=1, at = loc_2, font=2, line=1)
text(x=mean(1:nrow(test))
     , y=midline-.0025, "coverage probability (95% desired)", col="black", font=2)
text(x=mean(1:nrow(test))
     , y=topline-0.0025, "RMSE", col="black", font=2)
legend(x=6.5, y=0.02, col = "#DF8F44", legend = "pixel-level specifications", seg.len=0.65, inset = 0.005,  box.lty=0, cex=0.8, lty = 1, lwd = 4, bg="transparent")


####################################################################################
#### 3 panel spec chart
####################################################################################


library(tidyverse)
summary_full <- readRDS("summary_full2.rds")

full_summary <- summary_full %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes) & pixel.fe == FALSE)%>%
  group_by(std_a, std_v, std_p, pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county)%>%
  summarize(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  dplyr::select(Bias, everything())%>%
  ungroup()



full_summary <- full_summary %>%  
  dplyr::arrange(std_p, abs(Bias), by.group=TRUE)

library(data.table)
full_summary <- setDT(full_summary)[full_summary[, c(.I, NA), std_p]$V1][!.N]

par(oma=c(1,0,1,1))

labels <- list(#"Model:" = " ",
  "Unit of analysis:" = c("pixel", "grid", "property", "county"),
  "Fixed effects:" = c("grid FE", "property FE", "county FE", "treatment FE"),
  "Weights:" = c("unit area"),
  "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))


f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))

coverage <- full_summary$cover
coverage[is.na(coverage)] <- 0
c_print <- f(round(full_summary$cover, digits = 2)*100)
c_print[is.na(c_print)] <- " "

RMSE <- full_summary$RMSE
RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)

RMSE_print<- f(round(full_summary$RMSE, digits =4))
RMSE_print[is.na(RMSE_print)] <- " "

column_indic <- c(1:12, " ", 13:22, " ", 23:32)

test <- as.data.frame(subset(full_summary, select=-c(pixel.fe, cover, RMSE, std_p, std_v, std_a)))%>%
  filter(pixel=="TRUE"|pixel=="FALSE")%>%
  mutate_at(vars(Bias, q05, q95), as.numeric)

highlight_rows <- which(test[ , "property.fe"] == TRUE )

index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart

topline = -0.03
midline = topline-0.03
ylim <- c(midline-0.015,0.0475)
#Create the plot

schart(test,labels, ylim = ylim, index.ci=index.ci, ylab="Bias", highlight=highlight_rows,
       col.est=c("black","#00A1D5"),
       col.dot=c("black","lightgrey","red","#00A1D5"),
       bg.dot=c("black","lightgrey","white","#00A1D5"),
       #,band.ref=c(-.05, .04)
       axes = FALSE, n=12
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = 0, labels=TRUE)
abline(h=topline)
abline(h=midline)
abline(v=13, lty="dashed")
abline(v=26, lty="dashed")
lapply(1:length(RMSE), function(i) {
  rect(xleft=i-.3, ybottom=midline, xright=i+.3, ytop=midline+RMSE[i]*2.2, border=NA, col="#D55E00")
  #mtext(paste0(column_indic[i]), side=1, at = i, font=2, cex=.9)#, line=1, at=-1)
  text(x= i, y=midline+RMSE[i]*2.2+0.004, paste0(RMSE_print[i]), col="black", font=1, cex=.65)
  text(x= i, y=midline-0.01, paste0(c_print[i]), col="black", font=1, cex=.75 )
})
# text(x=5#mean(1:nrow(test))
#      , y=topline-0.0075, "RMSE", col="black", font=2)
mtext("RMSE", side=2, at = midline+0.015, font=2, las=1, line=.5)
# text(x=5#mean(1:nrow(test))
#      , y=midline-0.0075, "coverage probability (%)", col="black", font=2)
mtext("Coverage\nprobability (%)", side=2, at = midline-0.01, font=2, las=1, line=.5)
text(x=2
     , y=.0475, expression(paste(sigma[p],"=0.0")), col="black", font=2, cex = 1.1)
text(x=15.5
     , y=.0475, expression(paste(sigma[p],"=0.1")), col="black", font=2, cex = 1.1)
text(x=29
     , y=.0475, expression(paste(sigma[p],"=0.25")), col="black", font=2, cex = 1.1)
legend(x=-3, y=0.065, col = c("#00A1D5"), legend = c("specifications\nincorporating\nproperty"), inset = 0.005,  box.lty=0, cex=0.95
       ,  seg.len=0.25, lty = 1, horiz=TRUE, lwd = 4, bg="transparent")


####################################################################################
#### alternative parameterizations
####################################################################################
# activate which df to read in to change parameterization
#df_summary <- readRDS("summary_long_alt.rds") %>%
df_summary <- readRDS("summary_long_alt2.rds") %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes) & pixel.fe == FALSE)%>%
  group_by(pixel, grid, property, county, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county)%>%
  summarize(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  dplyr::select(Bias, everything())%>%
  ungroup()%>%
  dplyr::arrange(std_p, abs(Bias), by.group=TRUE)

#df_summary <- rbind(df_summary[11,], df_summary[4:10,], df_summary[1:3,])

par(oma=c(1,0,1,1))

labels <- list(#"Model:" = " ",
  "Unit of analysis:" = c("pixel", "grid", "property", "county"),
  "Fixed effects:" = c("grid FE", "property FE", "county FE", "treatment FE"),
  "Weights:" = c("unit area"),
  "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))

f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))
coverage <- f(round(df_summary$cover, digits=3)*100)
RMSE <- f(round(df_summary$RMSE, digits=4))
test <- as.data.frame(subset(df_summary, select=-c(cover, RMSE)))
highlight_rows <- which(test[ , "property.fe"] == TRUE )


index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart
ylim <- c(-0.04,0.03)

topline = -0.02

midline = topline-0.02-.0025

ylim <- c(midline-.02,0.04)
#bottomline = (min(ylim)+topline)/2
#Create the plot

schart(test,labels = labels, ylim = ylim, index.ci=index.ci, 
       col.est=c("black","#00A1D5"),
       col.dot=c("black","lightgrey","red","#00A1D5"),
       bg.dot=c("black","lightgrey","white","#00A1D5")
       , ylab="Bias", highlight=highlight_rows
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c( 0), labels=TRUE)
abline(h=topline)
abline(h=midline)
#abline(h=bottomline)
lapply(1:length(RMSE), function(i) {
  # text(x= i, y=min(ylim)+.002, paste0(i), col="black", font=2, cex = 1)
  #mtext(paste0(i), side=1, at = i, font=2, cex=.95)#, line=1, at=-1)
  text(x= i, y=midline+0.005, paste0(RMSE[i]), col="black", font=1, cex=.9)
  text(x= i, y=min(ylim)+0.0025, paste0(coverage[i]), col="black", font=1, cex=.9 )
})
# mtext('coverage\nprobability', side=2, las=1, at = loc_1, adj = 1, font=2, line=1)
# mtext('RMSE', side=2, las=1, at = loc_2, font=2, line=1)
text(x=mean(1:nrow(test))
     , y=midline-.005, "coverage probability (%)", col="black", font=2)
text(x=mean(1:nrow(test))
     , y=topline-.005, "RMSE", col="black", font=2)
legend(x=6.5, y=0.05, col = "#00A1D5", legend = "specifications\nincorporating\nproperty", seg.len=0.65, inset = 0.005,  box.lty=0, cex=0.8, lty = 1, lwd = 4, bg="transparent")

#########################################################################################################################################################
###### four panel spec chart
###########################################################################################################################################


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# summary_full.rds IS NOT CURRENT, USE summary_full2.rds and the 3 panel spec chart
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary_full <- readRDS("unbiased_dgp/summary_full.rds")

full_summary <- summary_full %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes) & pixel.fe == FALSE)%>%
  group_by(std_a, std_v, std_p, pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())%>%
  ungroup()



full_summary <- full_summary %>%  
  dplyr::arrange(std_p, abs(Bias), by.group=TRUE)

library(data.table)
full_summary <- setDT(full_summary)[full_summary[, c(.I, NA), std_p]$V1][!.N]

par(oma=c(1,0,1,1))

labels <- list(#"Model:" = " ",
  "Unit of analysis:" = c("pixel", "grid", "property", "county"),
  "Fixed effects:" = c("grid FE", "property FE", "county FE", "treatment FE"),
  "Weights:" = c("unit area"),
  "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))


f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))

coverage <- full_summary$cover
coverage[is.na(coverage)] <- 0
c_print <- f(round(full_summary$cover, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- full_summary$RMSE
RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)

RMSE_print<- f(round(full_summary$RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "

column_indic <- c(1:12, " ", 13:22, " ", 23:32, " ", 33:42, " ")

test <- as.data.frame(subset(full_summary, select=-c(pixel.fe, cover, RMSE, std_p, std_v, std_a)))%>%
  filter(pixel=="TRUE"|pixel=="FALSE")%>%
  mutate_at(vars(Bias, q05, q95), as.numeric)

highlight_rows <- which(test[ , "treatment.fe"] == FALSE )

index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart

topline = -0.03
midline = topline-0.03
ylim <- c(midline-0.015,0.0475)
#Create the plot

schart(test,labels, ylim = ylim, index.ci=index.ci, ylab="Bias", highlight=highlight_rows,
       col.est=c("black","#00A1D5"),
       col.dot=c("black","lightgrey","red","#00A1D5"),
       bg.dot=c("black","lightgrey","white","#00A1D5"),
       #,band.ref=c(-.05, .04)
       axes = FALSE, n=12
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = 0, labels=TRUE)
abline(h=topline)
abline(h=midline)
abline(v=13, lty="dashed")
abline(v=26, lty="dashed")
abline(v=39, lty="dashed")
lapply(1:length(RMSE), function(i) {
  rect(xleft=i-.3, ybottom=midline, xright=i+.3, ytop=midline+RMSE[i]*2.2, border=NA, col="#D55E00")
  #mtext(paste0(column_indic[i]), side=1, at = i, font=2, cex=.9)#, line=1, at=-1)
  #text(x= i, y=midline+RMSE[i]*1.5+0.004, paste0(RMSE_print[i]), col="black", font=1, cex=.65)
  text(x= i, y=midline-0.01, paste0(c_print[i]), col="black", font=1, cex=.75 )
})
# text(x=5#mean(1:nrow(test))
#      , y=topline-0.0075, "RMSE", col="black", font=2)
mtext("RMSE", side=2, at = midline+0.015, font=2, las=1, line=.5)
# text(x=5#mean(1:nrow(test))
#      , y=midline-0.0075, "coverage probability", col="black", font=2)
mtext("Coverage\nprobability", side=2, at = midline-0.01, font=2, las=1, line=.5)
text(x=1
     , y=.0475, expression(paste(sigma[p],"=0.0")), col="black", font=2, cex = 1.1)
text(x=15
     , y=.0475, expression(paste(sigma[p],"=0.1")), col="black", font=2, cex = 1.1)
text(x=28
     , y=.0475, expression(paste(sigma[p],"=0.2")), col="black", font=2, cex = 1.1)
text(x=41
     , y=.0475, expression(paste(sigma[p],"=0.3")), col="black", font=2, cex = 1.1)
legend(x=-3, y=0.06, col = c("#00A1D5"), legend = c("specifications\nincorporating\nspatial aggregation"), inset = 0.005,  box.lty=0, cex=0.95
       ,  seg.len=0.25, lty = 1, horiz=TRUE, lwd = 4, bg="transparent")

