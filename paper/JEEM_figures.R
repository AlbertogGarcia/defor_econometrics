# Load packages
library(knitr)
library(ggplot2)
library(kableExtra)
library(tidyverse)
library(Metrics)
library(reshape2)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggplot2)
library(scales)
library(data.table)

# schart function
source(here::here('paper', 'schart.R'))

# Color palette
palette <- list("white" = "#FAFAFA",
                "light_grey" = "#d9d9d9",
                "dark" = "#0c2230",
                "red" = "#ed195a",
                "blue" = "#1c86ee",
                "dark_blue" = "#00008b",
                "green" = "#00ab5b",
                "dark_green" = "#496F5D",
                "gold" = "#DAA520",
                "purple" = "#880ED4")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### set output directory
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_dir <- here::here("paper", "jeem_exhibits")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### set working directory
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setwd(here::here("paper"))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### TABLE 1: Lit table
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lit_table <- read.csv("JEEM_lit_table.csv")%>%#[,1:4]
  filter(Method != "matching")

kable(lit_table, format = "latex", row.names = FALSE,  booktabs = T, linesep = "",
      label = "table-lit",
      caption = "Example studies that measure conservation impacts by applying panel econometric methods to remotely sensed measures of deforestation.",
      col.names = c("Paper",
                    "Panel method",
                    "Unit of analysis",
                    "Unit FE level")) %>% 
  row_spec(0,bold=TRUE) %>% 
  kable_styling(font_size = 10, latex_options = c("HOLD_position"),
                #full_width = TRUE
                position = "left"
  )%>% 
  kableExtra::save_kable(paste0(out_dir, "/Table1.tex"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Table 2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# nobs in counterfactual landscape map
fig_nobs = 75^2

# Monte Carlo simulation parameters
nobs = 150^2
ppoints = 225
cpoints = 25
avg_parea = nobs/ppoints
avg_carea = nobs/cpoints
cellsize_med = 10
cellsize_large = 30

# property: 100 pixels | 100 x 30^2 = 90000 m^2 = 9 ha
p_ha = avg_parea * 30^2 / 10000
c_ha = avg_carea * 30^2 / 10000

# unit | grid cell resolution (m) | area (ha) | comparable use in literature

resolution_df <- data.frame(
  "unit" = c("Property", "County", "Large grid", "Small grid"),
  "avg_pix" = c(avg_parea, avg_carea, cellsize_large^2, cellsize_med^2),
  "structure" = c("Thiessen polygons", "Thiessen polygons", "Uniform square", "Uniform square")
) %>%
  mutate(area_ha = avg_pix * 30^2 / 10000) %>% 
  select(unit, structure, avg_pix, area_ha)

fig_nobs = 75^2
analysis_to_fig_ratio = nobs/fig_nobs

analysis_km2 <- nobs * 0.0009

kable(resolution_df, format = "latex", row.names = FALSE,  booktabs = T,
      label = "res",
      caption = "Spatial unit structure and size",
      col.names = NULL) %>%
  add_header_above(c("Spatial unit"=1, "Spatial structure"=1, "Avg. number of pixels"=1, "Area (hectares)"=1), escape = TRUE, line =TRUE, bold = TRUE)%>%
  kable_styling(font_size = 9.5, latex_options = c("HOLD_position"),
                position = "center")%>% 
  kableExtra::save_kable(paste0(out_dir, "/Table2.tex"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Fig 2 (Fig 1 is already saved)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results_long <- readRDS("results/results_aggregation.rds") 

df_summary <- results_long %>%
  filter(is.na(notes),
         weights == 0 ,
         (cox == 0 | HE.estimator == 1),
         (grid.fe == 0 | gridsize == max(gridsize, na.rm = T)),
         (property.fe == 0 | se_county == 0)
  ) %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, HE.estimator, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, HE.estimator, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())

df_summary <- rbind(df_summary[9,],
                    df_summary[5:8,], 
                    df_summary[1:3,],
                    df_summary[4,]
                    )


labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "Survival" = c("ATT-Cox estimator"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))

na_row <- setNames(data.frame(matrix(ncol = ncol(df_summary), nrow = 0)), colnames(df_summary))
na_row[1,] <- NA

select_results <- rbind(na_row, as.data.frame(df_summary), na_row)

coverage <- select_results$cover
c_print <- round(coverage, digits = 3)
c_print[is.na(c_print)] <- ""

RMSE <- select_results$RMSE
RMSE <- as.numeric(RMSE)
RMSE_print<- round(RMSE, digits = 4)
RMSE_print[is.na(RMSE_print)] <- ""

select_results <- as.data.frame(subset(select_results, select=-c(cover, RMSE)))

index.ci <- match(c("q05","q95"), names(select_results))

# One could also add information about model fit to this chart
topline = -0.007

midline = topline-0.004-.0022

ylim <- c(midline-.003,0.045)
#bottomline = (min(ylim)+topline)/2
#Create the plot

rmse_cex = 0.95
cov_cex = 1
leg_cex = 1

png(paste0(out_dir,"/Figure2.png"), width = 9, height = 11, units = "in", res = 300)

par(oma=c(1,0,1,1))

schart(select_results,labels, ylim = ylim, index.ci=index.ci, col.est = c(palette$dark, palette$red),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red),
       ylab="    Bias", highlight=c(2)
       , heights = c(2.5,1)
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c(-0.005, 0, 0.005, 0.02, 0.04), labels=TRUE)
abline(h=topline)
abline(h=midline)
#abline(h=bottomline)
lapply(1:(length(RMSE_print)), function(i) {
  mtext(if(i<=9){paste0(i)}, side=1, at = i + 1, font=2, cex=.95)#, line=1, at=-1)
  text(x= i, y=midline+0.0015, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=min(ylim)-0.0008, paste0(c_print[i]), col="black", font=1, cex=cov_cex )
})
segments(6.5, topline, x1 = 6.5, y1 = 0.041, lty = 2)
segments(9.5, topline, x1 = 9.5, y1 = 0.041, lty = 2)
text(x=mean(1:nrow(select_results))
     , y=midline - 0.0015, "Coverage probability", col="black", font=2)
text(x=mean(1:nrow(select_results))
     , y=topline-.0018, "RMSE", col="black", font=2)
text(x=4.5 , y=0.038, "Aggregated\nfixed effects", col=palette$dark, font=1)
text(x=7.9 , y=0.038, "Aggregated unit\nof analysis", col=palette$dark, font=1)
text(x=10.1 , y=0.0065, expression(widehat(ATT)["Cox"]), col=palette$dark, font=1, cex = 0.9)
text(x=2, y=0.0315, "Pixel-level\nTWFE", col=palette$red, font=1, cex = 0.9)
text(x=3, y=-0.0045, "pixel-level\nDID", col=palette$dark, font=1, cex = 0.9)
legend(x=8, y=0.047, col = "black", legend = "0.05 to 0.95 quantile\nof bias distribution", seg.len=0.65, inset = 0.005,  box.lty=0, cex=leg_cex, lty = 1, lwd = 4, bg="transparent")

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Fig 3 (Cox models)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_proportional <- readRDS("results/results_proportional.rds") %>% filter(cox == 1) %>% mutate(t_assumption = "proportional")
results_parallel <- readRDS("results/results_aggregation.rds") %>% filter(cox == 1)%>% mutate(t_assumption = "parallel")

df_summary <- rbind(results_parallel, results_proportional) %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  group_by(t_assumption, HE.estimator)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())

na_row <- setNames(data.frame(matrix(ncol = ncol(df_summary), nrow = 0)), colnames(df_summary))
na_row[1,] <- NA

select_results <- rbind(na_row, as.data.frame(df_summary)[1:2,], na_row, as.data.frame(df_summary)[3:4,], na_row)


coverage <- select_results$cover
c_print <- round(coverage, digits = 3)
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
RMSE <- as.numeric(RMSE)
RMSE_print<- round(RMSE, digits = 4)
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>%
  mutate(F_col = FALSE, F_col2 = FALSE)%>%
  select(c(Bias, F_col, F_col2, q05, q95))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("", "")

topline = -0.006
midline = topline-0.003
ylim <- c(midline-0.002,0.017)

png(paste0(out_dir,"/Figure_proportional.png"), width = 8, height = 5, units = "in", res = 300)

par(oma=c(1, 0,1,0.5))

schart(as.data.frame(schart_results), labels = labels, 
       ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="              Bias",
       col.est = c(palette$dark, palette$red), 
       bg.dot=c("white", "white", "white", "white"),
       col.dot=c("white", "white", "white", "white"),
       heights = c(100,1), cex = c(1.25,1.25))
Axis(side=2, at = c(-0.01, -0.005, 0, 0.005, 0.01), labels=TRUE)
segments(4, topline, x1 = 4, y1 = 0.019)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=topline-0.002, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=midline-0.0023, paste0(c_print[i]), col="black", font=1, cex=cov_cex)
})
text(x=mean(1:nrow(schart_results))
     , y=midline-.001, "Coverage probability", col="black", font=2)
text(x=mean(1:nrow(schart_results))
     , y=topline-.001, "RMSE", col="black", font=2)
text(x=2 , y=0.0115, expression(widehat(beta)["Cox"]), col=palette$dark, font=1, cex = 0.9)
text(x=3 , y=0.0053, expression(widehat(ATT)["Cox"]), col=palette$dark, font=1, cex = 0.9)
text(x=5 , y=0.0065, expression(widehat(beta)["Cox"]), col=palette$dark, font=1, cex = 0.9)
text(x=6 , y=0.0138, expression(widehat(ATT)["Cox"]), col=palette$dark, font=1, cex = 0.9)
text(x=1.75 , y=0.017, "Parallel trends holds", col=palette$dark, font=1)
text(x=5.15 , y=0.017, "Proportional trends holds", col=palette$dark, font=1)

dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Fig 4 (selection bias)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

specchart_long <- readRDS("results/results_selection.rds")


df_summary <- specchart_long %>%
  filter(is.na(notes),
         pixel.fe ==0,
         weights == 0 ,
         cox == 0,
         (grid.fe == 0 | gridsize == max(gridsize, na.rm = T)),
         (property.fe == 0 | se_county == 0)
  ) %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())

df_summary <- rbind(df_summary[4:7,],
                    df_summary[1:3,])


labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               #"Survival" = c("ATT-Cox estimator"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))

na_row <- setNames(data.frame(matrix(ncol = ncol(df_summary), nrow = 0)), colnames(df_summary))
na_row[1,] <- NA

select_results <- rbind(na_row, as.data.frame(df_summary), na_row)


coverage <- select_results$cover
c_print <- round(coverage, digits = 3)
c_print[is.na(c_print)] <- ""

RMSE <- select_results$RMSE
RMSE <- as.numeric(RMSE)
RMSE_print<- round(RMSE, digits = 4)
RMSE_print[is.na(RMSE_print)] <- ""

select_results <- as.data.frame(subset(select_results, select=-c(cover, RMSE)))

index.ci <- match(c("q05","q95"), names(select_results))

# One could also add information about model fit to this chart
topline = -0.007

midline = topline-0.004-.0025

ylim <- c(midline-.005,0.008)
#bottomline = (min(ylim)+topline)/2
#Create the plot

rmse_cex = 1
cov_cex = 1.1
leg_cex = 1.2

png(paste0(out_dir,"/Figure4.png"), width = 9, height = 8, units = "in", res = 300)

par(oma=c(1,0,1,1))

schart(select_results,labels, ylim = ylim, index.ci=index.ci, col.est = c(palette$dark, palette$dark_blue),
       bg.dot=c(palette$dark, "grey95", "white", palette$dark_blue),
       col.dot=c(palette$dark, "grey95", "white", palette$dark_blue),
       ylab="    Bias"#, highlight=c(6,7,8)
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       , heights = c(3,2)
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c(-0.005, 0, 0.005), labels=TRUE)
abline(h=topline)
abline(h=midline)
segments(5.5, topline, x1 = 5.5, y1 = 0.015, lty = 2)
lapply(1:(length(RMSE_print)), function(i) {
  text(x= i, y=midline+0.0015, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=min(ylim)+0.0008, paste0(c_print[i]), col="black", font=1, cex=cov_cex )
})
text(x=mean(1:nrow(select_results))
     , y=midline - 0.0015, "Coverage probability", col="black", font=2)
text(x=mean(1:nrow(select_results))
     , y=topline-.0018, "RMSE", col="black", font=2)
text(x=3.3 , y=0.0065, "Aggregated\nfixed effects", col=palette$dark, font=1)
text(x=7 , y=0.0065, "Aggregated unit\nof analysis", col=palette$dark, font=1)

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Fig 5 (DID property)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results_full <- readRDS("results/results_full.rds")

full_summary <- results_full %>%
  filter(treatment.fe == 1 & se_pixel == 1 & is.na(notes) & cox == 0) %>%
  mutate_at(vars(bias, cover), as.numeric) %>%
  group_by(std_p) %>%
  summarise(Bias = mean(bias),
            cover = mean(cover),
            RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95)) %>%
  mutate(s0 = std_p==0,
         s1 = std_p==0.1,
         s2 = std_p==0.2,
         s3 = std_p==0.3) %>% 
  arrange(desc(s0), desc(s1), desc(s2), desc(s3))

na_row <- setNames(data.frame(matrix(ncol = ncol(full_summary), nrow = 0)), colnames(full_summary))
na_row[1,] <- NA

select_results <- rbind(na_row, as.data.frame(full_summary), na_row)


coverage <- select_results$cover
c_print <- round(coverage, digits = 3)
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
RMSE <- as.numeric(RMSE)
RMSE_print<- round(RMSE, digits = 4)
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>%
  mutate(F_col = FALSE,
         F_col2 = FALSE)%>%
  select(c(Bias, F_col, F_col2, q05, q95))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("",
               "")


topline = -0.012
midline = topline-0.003
ylim <- c(midline-0.002,0.003)

png(paste0(out_dir,"/Figure5.png"), width = 10, height = 7, units = "in", res = 300)

par(oma=c(1,0,1,1))

schart(as.data.frame(schart_results), labels = labels, 
       ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="              Bias",
       #leftmargin = 16, 
       col.est = c(palette$dark, palette$red), 
       bg.dot=c("white", "white", "white", "white"),
       col.dot=c("white", "white", "white", "white"),
       heights = c(6,1), cex = c(1.25,1.25))

#mtext("Bias             ", side=2, at = -0.00, font=2, las=1, line=.5)
Axis(side=2, at = c(-0.01, -0.005, 0, 0.005, 0.01), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=topline-0.002, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=midline-0.002, paste0(c_print[i]), col="black", font=1, cex=cov_cex)
})
text(x=mean(1:nrow(schart_results))
     , y=midline-.001, "Coverage probability", col="black", font=2)
text(x=mean(1:nrow(schart_results))
     , y=topline-.001, "RMSE", col="black", font=2)
# mtext("RMSE", side=2, at = midline+0.003, font=2, las=1, line=.5)
# mtext("Coverage\nprobability", side=2, at = midline-0.003, font=2, las=1, line=.5)
mtext(paste0(0), side=1, at = 2, font=2, cex=1.2, line=3)
mtext(paste0(0.1), side=1, at = 3, font=2, cex=1.2, line=3)
mtext(paste0(0.2), side=1, at = 4, font=2, cex=1.2, line=3)
mtext(paste0(0.3), side=1, at = 5, font=2, cex=1.2, line=3)
mtext(bquote("Value of " ~ sigma[p]), side=1, at = 0, font=2, cex=1.25, line=3)

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Fig 6 (aggregation property)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_p = 10
grid_c = 30
grid_small = grid_p
grid_large = grid_c



full_summary <- results_full %>%
  mutate_at(vars(bias, cover, power), as.numeric)%>%
  filter(std_p==0.3 , 
         is.na(notes) , 
         gridsize %in% c(grid_p, grid_c) | grid.fe == 0,
         pixel.fe == 0 , 
         weights == 0 , 
         cox == 0 
  ) %>%
  group_by(pixel, grid, property, county, grid.fe, property.fe, county.fe, treatment.fe, se_pixel, se_grid, se_property, se_county, gridsize)%>%
  summarise(Bias = mean(bias),
            cover = mean(cover),
            power = mean(power),
            RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95))%>%
  mutate(grid.fe_large = ifelse(gridsize == grid_large, 1, 0),
         grid.fe_small = ifelse(gridsize == grid_small, 1, 0),
         grid_large = ifelse(grid.fe_large == 1 & grid == 1, 1, 0),
         grid_small = ifelse(grid.fe_small == 1 & grid == 1, 1, 0)
  )%>%
  mutate_at(vars(grid.fe_small, grid.fe_large), ~replace(., is.na(.), 0))

# Panel A - Pixel-level analyses
select_results_A <- full_summary %>% 
  filter(pixel==1) %>% 
  dplyr::arrange(Bias)
# dplyr::arrange(property.fe, county.fe, grid.fe, gridsize)

select_results_B <- full_summary %>% 
  filter(pixel==0) %>% 
  dplyr::arrange(Bias)
# dplyr::arrange(property, county, grid, gridsize)

na_row <- setNames(data.frame(matrix(ncol = ncol(select_results_A), nrow = 0)), colnames(select_results_A))
na_row[1,] <- NA

select_results <- rbind(na_row, select_results_A, na_row, select_results_B, na_row)

coverage <- select_results$cover
c_print <- round(coverage, digits = 4)
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
RMSE <- as.numeric(RMSE)
RMSE_print<- round(RMSE, digits =4)
RMSE_print[is.na(RMSE_print)] <- " "

gridsize_print <- select_results$gridsize
gridsize_print[is.na(gridsize_print)] <- " "

schart_results <- select_results %>% 
  select(c(Bias, q05, q95, pixel, grid_small, grid_large, county, property, treatment.fe, grid.fe_small, grid.fe_large, county.fe, property.fe
  ))%>%
  mutate_at(vars(pixel, county, property, grid_small, grid_large, treatment.fe, county.fe, property.fe, grid.fe_large, grid.fe_small), ~as.logical(as.numeric(.)))

index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("Unit of analysis:" = c("pixel", "small grid", "large grid", "county", "property"),
               "Fixed effects:" = c("treatment FE", "small grid FE", "large grid FE", "county FE", "property FE")
)

topline = -0.0125
midline = topline-0.005
ylim <- c(midline-0.004,0.0078)

png(paste0(out_dir,"/Figure6.png"), width = 10, height = 7, units = "in", res = 300)

par(oma=c(1,0,1,1))

schart(as.data.frame(schart_results), labels = labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="        Bias", 
       highlight = c(6, 11),  
       #  leftmargin = 16, 
       col.est = c(palette$dark, palette$red), heights = c(6,4), cex = c(1.25, 1.25),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red)
)
text(x=7
     , y=midline-.0015, "Coverage probability", col="black", font=2)
text(x=7
     , y=topline-.0015, "RMSE", col="black", font=2)
Axis(side=2, at = c(-0.01, 0.005, 0, 0.005), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0015, paste0(RMSE_print[i], " "), col="black", font=1, cex=rmse_cex)
  text(x= i, y=midline-0.004, paste0(c_print[i]), col="black", font=1, cex=cov_cex)
})
text(x=4 , y=0.0065, "Aggregated\nfixed effects", col=palette$dark, font=1)
text(x=9.5 , y=0.0065, "Aggregated unit\nof analysis", col=palette$dark, font=1)

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Fig 7 (weighting)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

summary_pweights <- readRDS("results/results_pweights.rds")

summary <- summary_pweights %>%
  mutate_at(vars(estimate, cover, p_ATT), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  filter(is.na(notes) & pixel.fe == FALSE & (treatment.fe == T | property == T) & (se_property == T | se_pixel == T))%>%
  group_by(pixel, property, property.fe, treatment.fe, weights, se_property, se_pixel)%>%
  summarise(RMSE = rmse(estimate, given_ATT),
            q05 = quantile(estimate, probs = .05),
            q95 = quantile(estimate, probs = .95),
            Estimate = mean(estimate),
            cover = mean(cover),
            ls_ATT = mean(ls_ATT),
            p_ATT = mean(p_ATT))%>%
  select(Estimate, everything())

ls_ATT <- round(mean(summary$ls_ATT), digits = 4)
p_ATT <- round(mean(summary$p_ATT), digits = 4)

df_summary <- summary 

png(paste0(out_dir,"/Figure7.png"), width = 8, height = 5, units = "in", res = 300)

par(oma=c(1,0,1,1))

labels <- list("Unit of analysis:" = c("pixel", "property"),
               "Fixed effects:" = c("property FE", "treatment FE"),
               "Weights:" = c("area weights"))
coverage <- round(df_summary$cover, digits=3)
RMSE <- round(df_summary$RMSE, digits=5)
select_results <- as.data.frame(subset(df_summary, select=-c(cover, RMSE)))

select_results <- subset(select_results, select=-c(ls_ATT, p_ATT, se_property, se_pixel))%>%
  distinct()

na_row <- setNames(data.frame(matrix(ncol = ncol(select_results), nrow = 0)), colnames(select_results))
na_row[1,] <- NA

select_results <- rbind(na_row, select_results, na_row)


#identifying quantile columns
index.ci <- match(c("q05","q95"), names(select_results))

# limits
ylim <- c(-0.018, -0.006)
#Create the plot

schart(select_results,labels, ylim = ylim,leftmargin = 12, heights = c(4, 3), index.ci=index.ci, col.est = c(palette$dark, palette$red),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red)
       , ylab="Coefficient estimate"
       , axes = F
) # make some room at the bottom
Axis(side=2, at = c(-0.006, -0.01, -0.014), labels=TRUE)
abline(h=mean(summary$p_ATT), col = palette$blue, lty = "dashed")
abline(h=mean(summary$ls_ATT), col = palette$dark_blue, lty = "dashed")
text(x=4.65 , y=mean(summary$p_ATT), "property-level\nATT", col=palette$blue, font=1)
text(x=4.65 , y=mean(summary$ls_ATT), "landscape\nATT", col=palette$dark_blue, font=1)

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### Fig 8 (observed multi1)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rate_landscape <- readRDS("results_multi/landscape.rds")

ggplot(data=rate_landscape, aes(x=year, y=defor, colour=Group))+
  geom_line(linewidth=1.5)+
  ylab("Deforestation rate")+
  xlab("Year")+
  scale_y_continuous(labels = scales::percent)+
  scale_color_manual(values=c(palette$dark, palette$blue, palette$red))+
  ggtitle("Annual deforestation across multiple groups")+
  theme_minimal(base_size = 14)+
  theme(#legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(.5, "cm"),
    legend.title.align = 0.5,
    text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggsave(filename = paste0(out_dir, "/Figure8.png"), width = 9, height = 5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### FIGURE 9
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

county_es <- readRDS("results_multi/county_long.rds") %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

pixel_es <- readRDS("results_multi/pixel_long.rds")%>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

color_scale = c("TWFE" = palette$red, "Gardner (2021)" = palette$dark, 
                "Callaway and Sant'Anna (2020)" = palette$blue, 
                "Borusyak, Jaravel, Spiess (2021)" = palette$dark_blue,  
                "Roth and Sant'Anna (2021)" = palette$purple, 
                "Sun and Abraham (2020)" = palette$gold, 
                "Truth" = palette$green)

# Get list of estimators
estimators = unique(county_es$estimator)

# Subset factor levels
levels = c("TWFE", 
           "Borusyak, Jaravel, Spiess (2021)", 
           "Callaway and Sant'Anna (2020)", 
           "Gardner (2021)",
           "Roth and Sant'Anna (2021)", 
           "Sun and Abraham (2020)",
           "Truth")
levels = levels[levels %in% estimators]

# Subset color scales
color_scale = color_scale[names(color_scale) %in% estimators]

out_leg = county_es %>%
  dplyr::mutate(
    ci_lower = q05,
    ci_upper = q95,
    estimator = factor(estimator, levels = levels)
  )

p_leg <- ggplot(out_leg, ggplot2::aes(x = term, y = estimate, color = estimator)) +
  geom_point(size = 2.6) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = color_scale)+
  labs(color = "Estimator") +
  guides(
    color = ggplot2::guide_legend(title.position = "top", nrow = 3)
  )

# Add CI and estimator level for plotting

out_pix = pixel_es %>%
  dplyr::mutate(
    ci_lower = q05,
    ci_upper = q95,
    estimator = factor(estimator, levels = levels)
  )

out_county = county_es %>%
  dplyr::mutate(
    ci_lower = q05,
    ci_upper = q95,
    estimator = factor(estimator, levels = levels)
  )%>% 
  filter(estimator != "Roth and Sant'Anna (2021)")


# position 
position = position_dodge(width = 0.5)

ylimmin = min(min(out_pix$ci_lower), min(out_county$ci_lower))
ylimmax = max(max(out_pix$ci_upper), max(out_county$ci_upper))
ylim = c(ylimmin, ylimmax)

p1 <- ggplot(out_pix, ggplot2::aes(x = term, y = estimate, color = estimator, ymin = ci_lower, ymax = ci_upper)) +
  geom_point(position = position, size = 2.6) +
  geom_errorbar(position = position) +
  geom_vline(xintercept = -0.5, linetype = "dashed", linewidth = 0.25) +
  labs(y = "Mean point estimate", x = "Event Time", color = "Estimator") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = color_scale) +
  guides(
    color = ggplot2::guide_legend(title.position = "top", nrow = 2)
  ) +
  ylim(ylimmin, ylimmax)+
  ggtitle(" Pixel as unit of analysis")+
  theme(plot.title = element_text(hjust = 0))

p2 <- ggplot(out_county, ggplot2::aes(x = term, y = estimate, color = estimator, ymin = ci_lower, ymax = ci_upper)) +
  geom_point(position = position, size = 2.6) +
  geom_errorbar(position = position) +
  geom_vline(xintercept = -0.5, linetype = "dashed", linewidth = 0.25) +
  labs(y = "Mean point estimate", x = "Event Time") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = color_scale) +
  ylim(ylimmin, ylimmax)+
  ggtitle("Aggregated unit of analysis")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0))

p2$labels$y <- " "

png(paste0(out_dir,"/Figure9.png"), width = 9, height = 6, units = "in", res = 300)

ggarrange(p1, p2, ncol=2, nrow=1, 
          legend.grob = get_legend(p_leg), 
          legend="bottom")

x_u = 0.535
y_u = 0.615
ht = 0.645
grid.rect(x= unit(x_u, "npc"), y = unit(y_u, "npc"), width = 0.9, height = ht, gp = gpar(lwd = 2, col = "black", fill = NA))
grid.rect(x= unit(0.5, "npc"), y = unit(y_u, "npc"), width = 0, height = ht, gp = gpar(lwd = 2, col = "black", fill = NA))

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### FIGURE 10
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
het_landscape <- readRDS("results_multi/het_landscape.rds")

ggplot(data=het_landscape, aes(x=year, y=defor, colour=Group))+
  geom_line(linewidth=1.5)+
  ylab("deforestation rate")+
  scale_y_continuous(labels = scales::percent)+
  scale_color_manual(values=c(palette$dark, palette$blue, palette$red))+
  ggtitle("Annual deforestation across heterogeneous groups")+
  theme_minimal(base_size = 14)+
  theme(#legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(.5, "cm"),
    legend.title.align = 0.5,
    text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
ggsave(filename = paste0(out_dir, "/Figure10.png"), width = 9, height = 5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### FIGURE 11
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
county_es <- readRDS("results_multi/county_long_hetTE.rds") %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

pixel_es <- readRDS("results_multi/pixel_long_hetTE.rds")%>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)


# create confidence intervals
out_pix = pixel_es %>%
  dplyr::mutate(
    ci_lower = q05,
    ci_upper = q95,
    estimator = factor(estimator, levels = levels)
  )

out_county = county_es %>%
  dplyr::mutate(
    ci_lower = q05,
    ci_upper = q95,
    estimator = factor(estimator, levels = levels)
  )%>% 
  filter(estimator != "Roth and Sant'Anna (2021)")


# position 
position = position_dodge(width = 0.5)

ylimmin = min(min(out_pix$ci_lower), min(out_county$ci_lower))
ylimmax = max(max(out_pix$ci_upper), max(out_county$ci_upper))
ylim = c(ylimmin, ylimmax)

p1 <- ggplot(out_pix, ggplot2::aes(x = term, y = estimate, color = estimator, ymin = ci_lower, ymax = ci_upper)) +
  geom_point(position = position, size = 2.6) +
  geom_errorbar(position = position) +
  geom_vline(xintercept = -0.5, linetype = "dashed", linewidth = 0.25) +
  labs(y = "Mean point estimate", x = "Event Time", color = "Estimator") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = color_scale) +
  guides(
    color = ggplot2::guide_legend(title.position = "top", nrow = 3)
  ) +
  ylim(ylimmin, ylimmax)+
  ggtitle(" Pixel as unit of analysis")+
  theme(plot.title = element_text(hjust = 0))

p2 <- ggplot(out_county , ggplot2::aes(x = term, y = estimate, color = estimator, ymin = ci_lower, ymax = ci_upper)) +
  geom_point(position = position, size = 2.6) +
  geom_errorbar(position = position) +
  geom_vline(xintercept = -0.5, linetype = "dashed", linewidth = 0.25) +
  labs(y = "Mean point estimate", x = "Event Time") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = color_scale) +
  ylim(ylimmin, ylimmax)+
  ggtitle("Aggregated unit of analysis")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0))

p2$labels$y <- " "

png(paste0(out_dir,"/Figure11.png"), width = 9, height = 6, units = "in", res = 300)

ggarrange(p1, p2, ncol=2, nrow=1, legend.grob = get_legend(p_leg), legend="bottom")

x_u = 0.535
y_u = 0.615
ht = 0.645
grid.rect(x= unit(x_u, "npc"), y = unit(y_u, "npc"), width = 0.9, height = ht, gp = gpar(lwd = 2, col = "black", fill = NA))
grid.rect(x= unit(0.5, "npc"), y = unit(y_u, "npc"), width = 0, height = ht, gp = gpar(lwd = 2, col = "black", fill = NA))

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Beginning Appendix
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Tile fig
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Setup the data
m <- matrix(c("Pixel 1:\nnever\ncleared","Pixel 2:\ncleared\nYear 2","Pixel 3:\ncleared\nYear 3","Pixel 4:\ncleared\nYear 1","Pixel 5:\ncleared\nYear 4","Pixel 6:\ncleared\nYear 2","Pixel 7:\ncleared\nYear 4","Pixel 8:\nnever\ncleared", "Pixel 9:\ncleared\nYear 3"), nrow=3, ncol=3, byrow = T)
df <- expand.grid(x=1:ncol(m),y=1:nrow(m))
df$val <- m[as.matrix(df[c('y','x')])]

ggplot(df, aes(x=x, y=y, label=val)) + 
  geom_tile(fill='transparent', colour = 'black') + 
  geom_text(size = 3) + 
  scale_y_reverse() +
  theme_classic() + 
  theme(axis.text  = element_blank(),
        panel.grid = element_blank(),
        axis.line  = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())+
  coord_fixed()
ggsave(filename = paste0(out_dir, "/A-tilefig.png"), width = 2.5, height = 2.5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### TWFE comp table
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

models_of_interest <- c("DID", "TWFE", "TWFE on dataset removing any pixel deforested pre-treatment", "ex-post difference in means")#, "final year difference in means")

twfe_comp <- readRDS("results/TWFE_comp.rds")%>%
  mutate(q25 = round(q25, digits = 5),
         q75 = round(q75, digits = 5),
         Bias = round(Bias, digits = 5),
         RMSE = round(RMSE, digits = 5),
         Model = ifelse(model == "TWFE on dataset dropping deforested pixels prior to treatment", "TWFE on dataset removing any pixel deforested pre-treatment", model),
         Model = ifelse(Model == "final period difference in means", "final year difference in means", Model))%>%
  unite("0.25 to 0.75 quantile", c(q25, q75), sep = " , ", remove = TRUE)%>%
  filter(Model %in% models_of_interest)%>%
  select(Model, Bias, everything(), - model)

kable(twfe_comp, format = "latex", row.names = FALSE,  booktabs = T,
      label = "twfe-comp",
      caption = "TWFE with pixel fixed effects is numerically equivalent to TWFE on dataset with all pixels deforested pre-treatment removed completely from the dataset"
      ,
      col.names = NULL
) %>%
  add_header_above(c("Model"=1, "Bias" = 1, "RMSE" = 1, "0.25 to 0.75 quantile" = 1), escape = TRUE, line =TRUE, bold = TRUE)%>%
  kable_styling(font_size = 10, latex_options = c("HOLD_position"),
                position = "center")%>% 
  kableExtra::save_kable(paste0(out_dir, "/A-twfe_comp.tex"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### complete summary fig
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results_full <- readRDS("results/results_full.rds")

full_summary <- results_full %>%
  filter(is.na(notes), 
         pixel.fe == 0, 
         (cox == 0 | HE.estimator ==1),
         (grid.fe == 0 | gridsize == max(gridsize, na.rm = T))
  ) %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, se_pixel, se_grid, se_property, se_county, HE.estimator), ~ as.logical(as.integer(.))
  )%>%
  group_by(std_a, std_v, std_p, pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, weights, HE.estimator, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())%>%
  ungroup()



full_summary <- full_summary %>%  
  dplyr::arrange(std_p, abs(Bias), by.group=TRUE)

full_summary <- setDT(full_summary)[full_summary[, c(.I, NA), std_p]$V1][!.N]

png(paste0(out_dir,"/A-full_summary.png"), width = 12, height = 12, units = "in", res = 300)

par(oma=c(1,0,1,1))

labels <- list(#"Model:" = " ",
  "Unit of analysis:" = c("pixel", "grid", "property", "county"),
  "Fixed effects:" = c("grid FE", "property FE", "county FE", "treatment FE"),
  # "Weights:" = c("unit area"),
  "Survival:" = c("ATT-Cox estimator"),
  "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))


coverage <- full_summary$cover
#coverage[is.na(coverage)] <- 0
c_print <- round(full_summary$cover, digits = 2)
c_print[is.na(c_print)] <- " "

RMSE <- full_summary$RMSE
RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)

RMSE_print<- round(full_summary$RMSE, digits =3)
RMSE_print[is.na(RMSE_print)] <- " "

select_results <- as.data.frame(subset(full_summary, select=-c(std_a, std_v, std_p, cover, RMSE, 
                                                               weights, pixel.fe)))

index.ci <- match(c("q05","q95"), names(select_results))

topline = -0.025

ylim <- c(topline,0.011)
#Create the plot

schart(select_results,labels, ylim = ylim, index.ci=index.ci, ylab="Bias", #highlight=highlight_rows,
       leftmargin = 15, 
       heights = c(5, 5), cex = c(1.75, 1.75),
       col.est = c(palette$dark, palette$red),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red),
       axes = FALSE#, n=8
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c(-.02, -.01, 0, 0.01), labels=TRUE)
#abline(h=topline)
#abline(h=midline)
abline(v=9, lty="dashed")
abline(v=18, lty="dashed")
abline(v=27, lty="dashed")
#lapply(1:length(RMSE), function(i) {
#rect(xleft=i-.3, ybottom=midline, xright=i+.3, ytop=midline+RMSE[i]*1.5, border=NA, col="#D55E00")
#mtext(paste0(column_indic[i]), side=1, at = i, font=2, cex=.9)#, line=1, at=-1)
#text(x= i, y=midline+RMSE[i]*1.5+0.004, paste0(RMSE_print[i]), col="black", font=1, cex=.65)
#  text(x= i, y=midline-0.01, paste0(c_print[i]), col="black", font=1, cex=.75 )
#})
# text(x=5#mean(1:nrow(test))
#      , y=topline-0.0075, "RMSE", col="black", font=2)
#mtext("RMSE", side=2, at = midline+0.015, font=2, las=1, line=.5)
# text(x=5#mean(1:nrow(test))
#      , y=midline-0.0075, "Coverage probability", col="black", font=2)
#mtext("Coverage\nprobability", side=2, at = midline-0.01, font=2, las=1, line=.5)
text(x=2
     , y=.01, expression(paste(sigma[p],"=0.0")), col="black", font=2, cex = 1.5)
text(x=11
     , y=.01, expression(paste(sigma[p],"=0.1")), col="black", font=2, cex = 1.5)
text(x=20
     , y=.01, expression(paste(sigma[p],"=0.2")), col="black", font=2, cex = 1.5)
text(x=29
     , y=.01, expression(paste(sigma[p],"=0.3")), col="black", font=2, cex = 1.5)
#legend(x=-3, y=0.06, col = c("#00A1D5"), legend = c("specifications\nincorporating\nspatial aggregation"), inset = 0.005,  box.lty=0, cex=0.95
#       ,  seg.len=0.25, lty = 1, horiz=TRUE, lwd = 4, bg="transparent")
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Alternate parameterizations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df_summary <- readRDS("results/results_alt.rds") %>%
  filter(is.na(notes),
         pixel.fe ==0,
         weights == 0 ,
         (cox == 0 | HE.estimator ==1),
         (grid.fe == 0 | gridsize == max(gridsize, na.rm = T)),
         (property.fe == 0 | se_county == 0)
  ) %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, HE.estimator, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, HE.estimator, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())

df_summary <- rbind(df_summary[5:8,],
                    df_summary[1:3,],
                    df_summary[4,])


labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "Survival" = c("ATT-Cox estimator"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))


select_results <- as.data.frame(df_summary)

cov_cex = 1
coverage <- select_results$cover
c_print <- round(coverage, digits = 3)
c_print[is.na(c_print)] <- ""

RMSE <- select_results$RMSE
RMSE <- as.numeric(RMSE)
RMSE_print<- round(RMSE, digits = 3)
RMSE_print[is.na(RMSE_print)] <- ""

select_results <- as.data.frame(subset(select_results, select=-c(cover, RMSE)))


highlight_rows <- which(select_results[ , "property.fe"] == TRUE )
index.ci <- match(c("q05","q95"), names(select_results))
#%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topline = -0.01
midline = topline-0.01-.0025
ylim <- c(midline-.01,0.035)

png(paste0(out_dir,"/A-alt_param1.png"), width = 9, height = 6, units = "in", res = 300)

par(oma=c(1,0,1,1))

schart(select_results,labels = labels, ylim = ylim, index.ci=index.ci, 
       col.est = c(palette$dark, palette$red),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red),
       ylab="Bias", highlight=highlight_rows
       ,band.ref=c(-.05, .04)
       , axes = FALSE
) # make some room at the bottom
Axis(side=2, at = c( 0, 0.02), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:(length(RMSE_print)), function(i) {
  text(x= i, y=midline+0.0024, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=min(ylim)+0.001, paste0(c_print[i]), col="black", font=1, cex=cov_cex )
})
text(x=mean(1:nrow(select_results))
     , y=midline-.003, "Coverage probability", col="black", font=2)
text(x=mean(1:nrow(select_results))
     , y=topline-.003, "RMSE", col="black", font=2)
legend(x=4, y=0.0365, col = palette$red, legend = "property-level specifications", seg.len=0.65, inset = 0.005,  box.lty=0, cex=1, lty = 1, lwd = 4, bg="transparent")
dev.off()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
######### creating chart for other alternate parameterization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

df_summary <- readRDS("results/results_alt2.rds") %>%
  filter(is.na(notes),
         pixel.fe ==0,
         weights == 0 ,
         (cox == 0 | HE.estimator ==1),
         (grid.fe == 0 | gridsize == max(gridsize, na.rm = T)),
         (property.fe == 0 | se_county == 0)
  ) %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, HE.estimator, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, HE.estimator, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())

df_summary <- rbind(df_summary[5:8,],
                    df_summary[1:3,],
                    df_summary[4,])


labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "Survival" = c("ATT-Cox estimator"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))


select_results <- as.data.frame(df_summary)

coverage <- select_results$cover
c_print <- round(coverage, digits = 3)
c_print[is.na(c_print)] <- ""

RMSE <- select_results$RMSE
RMSE <- as.numeric(RMSE)
RMSE_print<- round(RMSE, digits = 3)
RMSE_print[is.na(RMSE_print)] <- ""

select_results <- as.data.frame(subset(select_results, select=-c(cover, RMSE)))


highlight_rows <- which(select_results[ , "property.fe"] == TRUE )
index.ci <- match(c("q05","q95"), names(select_results))
#%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topline = -0.037
midline = topline-0.01-.0028
ylim <- c(midline-.01,0.01)

png(paste0(out_dir,"/A-alt_param2.png"), width = 9, height = 6, units = "in", res = 300)

par(oma=c(1,0,1,1))

schart(select_results,labels = labels, ylim = ylim, index.ci=index.ci, 
       col.est = c(palette$dark, palette$red),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red), 
       ylab="Bias", highlight=highlight_rows
       ,band.ref=c(-.05, .04)
       , axes = FALSE
) # make some room at the bottom
Axis(side=2, at = c( -0.02, 0), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:(length(RMSE_print)), function(i) {
  text(x= i, y=midline+0.0024, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=min(ylim)+0.001, paste0(c_print[i]), col="black", font=1, cex=cov_cex )
})
text(x=mean(1:nrow(select_results))
     , y=midline-.003, "Coverage probability", col="black", font=2)
text(x=mean(1:nrow(select_results))
     , y=topline-.003, "RMSE", col="black", font=2)
dev.off()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Defor rates
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outcome <- readRDS("results/outcomes.rds")
suppressWarnings(cbias <- melt(outcome, value.name = "bias"))

caption <- paste0("Outcome 1. Mean bias: ", round(mean(outcome$outcome1), digits = 4),"; RMSE: ",round(rmse(0, outcome$outcome1), digits = 5), "\n", "Outcome 2. Mean bias: ", round(mean(outcome$outcome2), digits = 4),"; RMSE: ",round(rmse(0, outcome$outcome2), digits = 5), "\n", "Outcome 3. Mean bias: ", round(mean(outcome$outcome3), digits = 4),"; RMSE: ",round(rmse(0, outcome$outcome3), digits = 5))

ggplot(data = cbias, aes(x = bias, fill=variable)) +
  geom_density(alpha = .6) +
  guides(fill=guide_legend(title=NULL))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x= "Bias", caption = caption) +
  scale_x_continuous(breaks = c(-0.01, -0.005, 0, 0.005))+
  theme_bw(base_size = 14)+
  theme(plot.caption = element_text(hjust = 0.5))+
  scale_fill_manual(labels=c("Outcome 1", "Outcome 2", "Outcome 3"), values = c(palette$dark, palette$red, palette$blue)
  )
ggsave(filename = paste0(out_dir, "/A-outcomes.png"), width = 7, height = 5)
