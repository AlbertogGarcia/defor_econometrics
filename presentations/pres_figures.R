
library(patchwork)
library(tidyverse)
library(ggplot2)
library(clubSandwich)
library(matrixStats)
library(ggplot2)
library(plm)
library(Metrics)
library(DataCombine)
library(dplyr)
library(DeclareDesign)
library(data.table)
source(here::here("unbiased_dgp", "schart.R"))
library(ggpubr)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helper functions --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_separated <- function(figure, filename){
  legend <- get_legend(figure)
  legend <- as_ggplot(legend)
  legend <- legend + 
    theme(plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
          panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
          plot.margin=unit(c(0,0.8,0,0),"in"))
  figure <- figure + theme(legend.position = "none")
  ggsave(paste0(out_dir, filename,".png"), figure, width = 3, height = 3, units = "in")
  ggsave(paste0(out_dir, filename, "_lgd", ".png"), legend, width = 3, height = 1.5, units = "in")
}

format_fig <- function(figure){
  figure <- figure +
    theme(panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
          plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.background=element_rect(fill = "#FAFAFA"),
          legend.key=element_rect(fill = "#FAFAFA"),
          legend.text=element_text(size=12))
  return(figure)
}

f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set parameters --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

palette <- list("white" = "#FAFAFA",
                "light_grey" = "#d9d9d9",
                "dark" = "#0c2230",
                "red" = "#ed195a",
                "blue" = "#1c86ee",
                "green" = "#7CAE7A",
                "dark_green" = "#496F5D",
                "gold" = "#DAA520")

results_dir <- paste0(getwd()[1], '/unbiased_dgp/')
out_dir <- paste0(getwd()[1], '/presentations/figs/')


base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

std_a = 0.1
std_v = 0.25
years = 3
nobs = 1000
n = 200

cellsize = 10
ppoints = 50
std_p = 0
cpoints = 20


std_avp = (std_a^2+std_v^2+std_p^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)


countyscape = full_landscape(nobs, cellsize, ppoints, cpoints)
pixloc_df = countyscape$pixloc_df
control_area = countyscape$control_area
intervention_area = countyscape$intervention_area
intervention_area_merge = intervention_area %>% st_union()
p_bounds = countyscape$p_bounds
c_bounds = countyscape$c_bounds


ATT <- pnorm(b0+b1+b2_1+b3, 0, (std_a^2+std_v^2 +std_p^2)^.5) - pnorm(b0+b1+b2_1, 0, (std_a^2+std_v^2 + std_p^2)^.5)

pixloc <- pixloc_df

Nobs <- length(pixloc$treat)  
panels <- fabricate(
  pixels = add_level(N = Nobs, a_i = rnorm(N, 0, std_a), treat = pixloc$treat),
  year = add_level(N = (years*2), nest = FALSE),
  obs = cross_levels(
    by = join(pixels, year),
    post = ifelse(year > years, 1, 0),
    v_it = rnorm(N, 0, std_v),
    ystar = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + b3*treat*post + a_i + v_it,
    ystar_counterfactual = b0 + b1*treat + b2_0*post*(1-treat) + b2_1*post*treat + a_i + v_it
  )
)


#generate random 
error_table <- data.frame(property = as.character(unique(pixloc$property)), p_err = rnorm(length(unique(pixloc$property)), 0, sd = std_p))

panels$pixels <- gsub("(?<![0-9])0+", "", panels$pixels, perl = TRUE)

panels <- panels %>%
  inner_join(pixloc, by = c("pixels", "treat")) %>%
  inner_join(error_table, by = "property") %>%
  mutate(ystar = ystar + p_err, 
         y = (ystar > 0)*1, 
         ystar_counterfactual = ystar_counterfactual+p_err,
         y_counterfactual = (ystar_counterfactual > 0)*1) 

panels_counterfactual <- panels

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
  mutate(indic = year - defor_year,
         defor = ifelse(indic > 0, 1, y),
         y_it = ifelse(indic > 0, NA, y) )



data_df <- panels %>%
  mutate(y_it = as.factor(replace_na(y_it, -1)),
         y_it = y_it %>% recode("-1" = "Previously deforested", 
                                "0" = "Not deforested", 
                                "1" = "Deforested"),
         treat = as.factor(treat),
         treat = treat %>% recode("0" = "Stable forest - not treated",
                                  "1" = "Stable forest - treated")) %>% 
  st_as_sf() %>%
  dplyr::select(pixels, year, treat, defor, y_it) 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Landscape plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fills = c("Previously deforested" = palette$light_grey,
          "Not deforested" = palette$dark, 
          "Deforested" = palette$red,
          "Stable forest - treated" = palette$dark,
          "Stable forest - not treated" = palette$blue)


# Pre-deforestation landscape
plot_df <- data_df %>%
  filter(year == 2) %>% 
  mutate(y_it = y_it %>% recode("Deforested" = "Not deforested"))

initial_forest <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = y_it), color = "white", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values= fills)
initial_forest <- format_fig(initial_forest)

initial_forest %>% save_separated("initial_forest")

# Period 2-6 deforestation
for (i in seq(2,6)){
  plot_df <- data_df %>%
    filter(year == i)
  defor_plot <- ggplot() + 
    geom_sf(data = plot_df, aes(fill = y_it), color = "white", shape = 22, alpha = 1, size = 2.8) +
    scale_fill_manual(values= fills)
  defor_plot %>% 
    format_fig %>% 
    save_separated(paste0("defor_", as.character(i))) 
}


# Visualizing exposure to treatment
plot_df <- data_df %>% 
  filter(year==2) %>% 
  mutate(plot_var = ifelse(y_it %in% c("Deforested", "Previously deforested"), as.character(y_it), as.character(treat)))
plot_df$plot_var <- factor(plot_df$plot_var, levels = c("Deforested", "Previously deforested", "Stable forest - treated", "Stable forest - not treated"))


intervention_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "white", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values = fills)
# geom_sf(data = c_bounds, aes(color = "county boundaries"), size = 1.5, fill = "NA", color = palette$white)
intervention_plot %>% 
  format_fig() %>% 
  save_separated(paste0("intervention")) 


# Treatment-level analysis
treat_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = intervention_area_merge, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
treat_plot %>% 
  format_fig() %>% 
  save_separated(paste0("treatment")) 


# County-level analysis
county_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = c_bounds, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
county_plot %>% 
  format_fig() %>% 
  save_separated(paste0("county")) 


# Grid-level analysis
grid_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values = fills)  + 
  geom_vline(xintercept = c(0,8,16,24,32), color="white", size=1.25) + 
  geom_hline(yintercept = c(0,8,16,24,32), color="white", size=1.25) +
  xlab("") +
  ylab("")


grid_plot <- grid_plot %>% 
  format_fig() %>% 
  save_separated(paste0("grid")) 


# Property-level analysis
prop_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = p_bounds, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
prop_plot %>% 
  format_fig() %>% 
  save_separated(paste0("property")) 





# # Add property boundaries
# prop_plot <- ggplot() + 
#   geom_sf(data = plot_df, aes(fill = plot_var, color = plot_var), shape = 22, alpha = 1, size = 2.8) +
#   scale_color_manual(values = fills) +
#   scale_fill_manual(values = fills) +
#   geom_sf(data = c_bounds, aes(color = "county boundaries"), size = 1.5, fill = "NA", color = palette$white) +
#   geom_sf(data = p_bounds, aes(color = "property boundaries"), fill = "NA", color = palette$gold)
# prop_plot %>% 
#   format_fig() %>% 
#   save_separated(paste0("property")) 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# TWFE vs DID plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twfe_df <- read.csv(paste0(results_dir, "TWFE_long2.csv")) %>% 
  select(-X)


ggplot() +
  xlim(c(-0.01, 0.05)) +
  ylim(c(0, 251)) +
  theme_bw(base_size = 10) +
  annotate("text", x = 0.04, y = 210, label = "Treated baseline = 0.03\nUntreated baseline = 0.02\nATT = -0.01\nTrend = -0.005", size = 3) +
  geom_vline(xintercept = 0, color = palette$red, linetype = "longdash", size = 0.5) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  xlab("Bias") +
  ylab("Frequency")
ggsave(paste0(out_dir, "twfe_0.png"), width = 5, height = 2.6, units = "in")

plot_df <- twfe_df %>% 
  filter(parameterization==2,
         model == "DID")
ggplot() +
  geom_density(data = plot_df, aes(x=bias, fill = model)) +
  xlim(c(-0.01, 0.05)) +
  ylim(c(0, 251)) +
  scale_fill_manual(values = c("DID"=palette$dark, "TWFE"=palette$blue)) +
  theme_bw(base_size = 10) +
  annotate("text", x = -0.0025, y = 25, color = palette$white, label = "DID", size = 3) +
  annotate("text", x = 0.04, y = 210, label = "Treated baseline = 0.03\nUntreated baseline = 0.02\nATT = -0.01\nTrend = -0.005", size = 3) +
  geom_vline(xintercept = 0, color = palette$red, linetype = "longdash", size = 0.5) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  xlab("Bias") +
  ylab("Frequency")
ggsave(paste0(out_dir, "twfe_1.png"), width = 5, height = 2.6, units = "in")


plot_df <- twfe_df %>% 
  filter(parameterization==2)
ggplot() +
  geom_density(data = plot_df, aes(x=bias, fill = model)) +
  xlim(c(-0.01, 0.05)) +
  ylim(c(0, 251)) +
  scale_fill_manual(values = c("DID"=palette$dark, "TWFE"=palette$blue)) +
  theme_bw(base_size = 10) +
  annotate("text", x = -0.0025, y = 25, color = palette$white, label = "DID", size = 3) +
  annotate("text", x = 0.0095, y = 25, color = palette$dark, label = "TWFE", size = 3) +
  annotate("text", x = 0.04, y = 210, label = "Treated baseline = 0.03\nUntreated baseline = 0.02\nATT = -0.01\nTrend = -0.005", size = 3) +
  geom_vline(xintercept = 0, color = palette$red, linetype = "longdash", size = 0.5) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  xlab("Bias") +
  ylab("Frequency")
ggsave(paste0(out_dir, "twfe_2.png"), width = 5, height = 2.6, units = "in")


plot_df <- twfe_df %>% 
  filter(parameterization==3)
ggplot() +
  geom_density(data = plot_df, aes(x=bias, fill = model)) +
  xlim(c(-0.01, 0.05)) +
  ylim(c(0, 251)) +
  scale_fill_manual(values = c("DID"=palette$dark, "TWFE"=palette$blue)) +
  theme_bw(base_size = 10) +
  annotate("text", x = -0.0025, y = 25, color = palette$white, label = "DID", size = 3) +
  annotate("text", x = 0.02, y = 25, color = palette$dark, label = "TWFE", size = 3) +
  annotate("text", x = 0.04, y = 210, label = "Treated baseline = 0.04\nUntreated baseline = 0.02\nATT = -0.01\nTrend = -0.005", size = 3) +
  geom_vline(xintercept = 0, color = palette$red, linetype = "longdash", size = 0.5) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  xlab("Bias") +
  ylab("Frequency")
ggsave(paste0(out_dir, "twfe_3.png"), width = 5, height = 2.6, units = "in")



plot_df <- twfe_df %>% 
  filter(parameterization==4)
ggplot() +
  geom_density(data = plot_df, aes(x=bias, fill = model)) +
  xlim(c(-0.01, 0.05)) +
  ylim(c(0, 251)) +
  scale_fill_manual(values = c("DID"=palette$dark, "TWFE"=palette$blue)) +
  theme_bw(base_size = 10) +
  annotate("text", x = -0.0025, y = 25, color = palette$white, label = "DID", size = 3) +
  annotate("text", x = 0.03, y = 25, color = palette$dark, label = "TWFE", size = 3) +
  annotate("text", x = 0.04, y = 210, label = "Treated baseline = 0.05\nUntreated baseline = 0.02\nATT = -0.01\nTrend = -0.005", size = 3) +
  geom_vline(xintercept = 0, color = palette$red, linetype = "longdash", size = 0.5) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  xlab("Bias") +
  ylab("Frequency")
ggsave(paste0(out_dir, "twfe_4.png"), width = 5, height = 2.6, units = "in")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spec charts: Emphasize aggregation as solution -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_summary <- read.csv(paste0(results_dir, "full_summary.csv")) %>% 
  select(-X)

## Option 1 - aggregate FEs
select_results <- full_summary %>% 
  filter(sigma_p==0,
         !((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         pixel==T) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel.fe), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))

na_row <- c("mean_bias" = NA,
            "q05" = NA,
            "q95" = NA)

select_results <- rbind(na_row, select_results, na_row)


coverage <- select_results$cover
# coverage[is.na(coverage)] <- 0
c_print <- f(round(coverage, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
# RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>% 
  # select(c(mean_bias, q05, q95, pixel, county, grid, property, treatment.fe, county.fe, grid.fe, property.fe))
  # select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, pixel.fe, treatment.fe, county.fe, grid.fe, property.fe))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("Pixel FE", "Treatment FE", "County FE", "Grid FE", "Property FE")


topline = -0.025
midline = topline-0.015
ylim <- c(midline-0.012,0.04)


png(paste0(out_dir,"spec_did_fe.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, labels = labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = 2, leftmargin = 20, col.est = c(palette$dark, palette$red), heights = c(3,1), cex = c(1.75,1.75))

mtext("Bias", side=2, at = 0.038, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.02, 0, 0.02), labels=TRUE, cex = 3)
abline(h=0, lty=2, lwd = 2)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=2)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=2)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5, cex = 2)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5, cex = 2)
dev.off()





## Option 2 - aggregate unit of observation
select_results <- full_summary %>% 
  filter(sigma_p==0,
         !((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         (pixel==F | pixel.fe==T)) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel.fe), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))

na_row <- c("mean_bias" = NA,
            "q05" = NA,
            "q95" = NA)

select_results <- rbind(na_row, select_results, na_row)


coverage <- select_results$cover
# coverage[is.na(coverage)] <- 0
c_print <- f(round(coverage, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
# RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>% 
  # select(c(mean_bias, q05, q95, pixel, county, grid, property, treatment.fe, county.fe, grid.fe, property.fe))
  # select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, pixel.fe, county.fe, grid.fe, property.fe))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("Pixel", "Grid", "County", "Property")


topline = -0.025
midline = topline-0.015
ylim <- c(midline-0.012,0.04)


png(paste0(out_dir,"spec_did_agg.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, labels = labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = 2, leftmargin = 20, col.est = c(palette$dark, palette$red), heights = c(3,1), cex = c(1.75,1.75))

mtext("Bias", side=2, at = 0.038, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.02, 0, 0.02), labels=TRUE, cex = 3)
abline(h=0, lty=2, lwd = 2)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=2)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=2)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5, cex = 2)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5, cex = 2)
dev.off()





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spec charts: Problem of property-level disturbances ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
select_results <- full_summary %>% 
  filter(treatment.fe==1,
         !((se_grid==1 | se_property==1 | se_county==1))) %>% 
  mutate(s0 = sigma_p==0,
         s1 = sigma_p==0.1,
         s2 = sigma_p==0.2,
         s3 = sigma_p==0.3) %>% 
  arrange(desc(s0), desc(s1), desc(s2), desc(s3))

na_row <- c("mean_bias" = NA,
            "q05" = NA,
            "q95" = NA)

select_results <- rbind(na_row, select_results, na_row)


coverage <- select_results$cover
# coverage[is.na(coverage)] <- 0
c_print <- f(round(coverage, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
# RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>% 
  # select(c(mean_bias, q05, q95, pixel, county, grid, property, treatment.fe, county.fe, grid.fe, property.fe))
  # select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, s0, s1, s2, s3))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("Value of sigma" = c("0", "0.1", "0.2", "0.3"))


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"spec_prop_did.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, labels = labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       leftmargin = 15, col.est = c(palette$dark, palette$red), heights = c(3,1), cex = c(1.75,1.75))

mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE, cex = 3)
abline(h=0, lty=2, lwd = 2)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=2)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=2)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5, cex = 2)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5, cex = 2)
dev.off()





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spec charts: structuring model around unit of analysis ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Panel A - Pixel-level analyses
select_results <- full_summary %>% 
  filter(sigma_p==0.3,
         !((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         !pixel.fe==T,
         pixel==T) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))

na_row <- c("mean_bias" = NA,
            "q05" = NA,
            "q95" = NA)

select_results <- rbind(na_row, select_results, na_row)


coverage <- select_results$cover
# coverage[is.na(coverage)] <- 0
c_print <- f(round(coverage, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
# RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>% 
  # select(c(mean_bias, q05, q95, pixel, county, grid, property, treatment.fe, county.fe, grid.fe, property.fe))
# select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, treatment.fe, county.fe, grid.fe, property.fe))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("Treatment FE", "County FE", "Grid FE", "Property FE")


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"spec_prop_fe.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, labels = labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = 5, leftmargin = 20, col.est = c(palette$dark, palette$red), heights = c(3,1), cex = c(1.75,1.75))

mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE, cex = 3)
abline(h=0, lty=2, lwd = 2)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=2)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=2)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5, cex = 2)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5, cex = 2)
dev.off()



# Panel B - Aggregate units of analysis
select_results <- full_summary %>% 
  filter(sigma_p==0.3,
         !((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         !pixel.fe==T,
         pixel==F) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))

na_row <- c("mean_bias" = NA,
            "q05" = NA,
            "q95" = NA)

select_results <- rbind(na_row, select_results, na_row)


coverage <- select_results$cover
# coverage[is.na(coverage)] <- 0
c_print <- f(round(coverage, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
# RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>% 
  # select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, county, grid, property))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("Grid", "County", "Property")


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"spec_prop_agg.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)

schart(schart_results, labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = 4, leftmargin = 25, col.est = c(palette$dark, palette$red), heights = c(3,1), cex = c(1.75, 1.75))

mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE, cex = 3)
abline(h=0, lty=2, lwd = 2)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=2)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=2)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5, cex = 2)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5, cex = 2)
dev.off()









# # Panel B - Aggregate units of analysis
# fe_results <- full_summary %>% 
#   filter(sigma_p==0.3,
#          !((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
#          !pixel.fe==T) %>% 
#   mutate(cover_dif = (cover - 0.95)*100) %>% 
#   arrange(desc(pixel), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))
# 
# 
# coverage <- fe_results$cover
# coverage[is.na(coverage)] <- 0
# c_print <- f(round(fe_results$cover, digits = 2))
# c_print[is.na(c_print)] <- " "
# 
# RMSE <- fe_results$RMSE
# RMSE[is.na(RMSE)] <- 0
# RMSE <- as.numeric(RMSE)
# RMSE_print<- f(round(fe_results$RMSE, digits =3))
# RMSE_print[is.na(RMSE_print)] <- " "
# 
# 
# schart_results <- fe_results %>% 
#   select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
#   select(c(mean_bias, q05, q95, pixel, county, grid, property, treatment.fe, county.fe, grid.fe, property.fe))
# index.ci <- match(c("q05","q95"), names(schart_results))
# 
# labels <- list(
#   "Unit of analysis:" = c("pixel", "county", "grid", "property"),
#   "Fixed effects:" = c("treatment FE", "grid FE", "property FE", "county FE"))
# 
# 
# topline = -0.03
# midline = topline-0.015
# ylim <- c(midline-0.012,0.02)
# 
# 
# # png(paste0(out_dir,"prop_spec.png"), width = 9, height = 7, units = "in", res = 150)
# par(bg = palette$white)
# 
# schart(schart_results, labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
#        highlight = c(4, 7), leftmargin = 12, col.est = c(palette$dark, palette$red))
# 
# mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5)
# 
# Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE)
# abline(h=topline)
# abline(h=midline)
# lapply(1:length(RMSE), function(i) {
#   text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=1)
#   text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=1)
# })
# 
# mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5)
# 
# mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5)
# dev.off()
# 


