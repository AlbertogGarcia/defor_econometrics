
library(patchwork)
library(tidyverse)

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


source(here::here("unbiased_dgp", "schart.R"))

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

nobs <- 1000

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


# Overlaying county boundaries
plot_df <- data_df %>%
  filter(year == 2)
county_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = y_it), color = "white", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = c_bounds, aes(color = "county boundaries"), size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
county_plot %>% 
  format_fig() %>% 
  save_separated(paste0("county")) 


# Visualizing exposure to treatment
plot_df <- data_df %>% 
  filter(year==2) %>% 
  mutate(plot_var = ifelse(y_it %in% c("Deforested", "Previously deforested"), as.character(y_it), as.character(treat)))


plot_df$plot_var <- factor(plot_df$plot_var, levels = c("Deforested", "Previously deforested", "Stable forest - treated", "Stable forest - not treated"))
intervention_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "white", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values = fills) +
  geom_sf(data = c_bounds, aes(color = "county boundaries"), size = 1.5, fill = "NA", color = palette$white)
intervention_plot %>% 
  format_fig() %>% 
  save_separated(paste0("intervention")) 


# Add property boundaries
prop_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var, color = plot_var), shape = 22, alpha = 1, size = 2.8) +
  scale_color_manual(values = fills) +
  scale_fill_manual(values = fills) +
  geom_sf(data = c_bounds, aes(color = "county boundaries"), size = 1.5, fill = "NA", color = palette$white) +
  geom_sf(data = p_bounds, aes(color = "property boundaries"), fill = "NA", color = palette$gold)
prop_plot %>% 
  format_fig() %>% 
  save_separated(paste0("property")) 



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





## Option 2 - aggregate unit of observation






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spec charts: Problem of property-level disturbances ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fe_results <- full_summary %>% 
  filter(!((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         !pixel.fe==T) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))


coverage <- fe_results$cover
coverage[is.na(coverage)] <- 0
c_print <- f(round(fe_results$cover, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- fe_results$RMSE
RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(fe_results$RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "


schart_results <- fe_results %>% 
  select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, pixel, county, grid, property, treatment.fe, county.fe, grid.fe, property.fe))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list(
  "Unit of analysis:" = c("pixel", "county", "grid", "property"),
  "Fixed effects:" = c("treatment FE", "grid FE", "property FE", "county FE"))


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"prop_spec.png"), width = 9, height = 7, units = "in", res = 150)
par(bg = palette$white)

schart(schart_results, labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = c(4, 7), leftmargin = 12, col.est = c(palette$dark, palette$red))

mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5)

Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=1)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=1)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5)
dev.off()





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spec charts: structuring model around unit of analysis ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Panel A - Pixel-level analyses
fe_results <- full_summary %>% 
  filter(sigma_p==0.3,
         !((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         !pixel.fe==T,
         pixel==T) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))


coverage <- fe_results$cover
coverage[is.na(coverage)] <- 0
c_print <- f(round(fe_results$cover, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- fe_results$RMSE
RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(fe_results$RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "


schart_results <- fe_results %>% 
  # select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, treatment.fe, county.fe, grid.fe, property.fe))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list(
  "Fixed effects:" = c("treatment FE", "grid FE", "property FE", "county FE"))


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"prop_spec.png"), width = 9, height = 7, units = "in", res = 150)
par(bg = palette$white)

schart(schart_results, labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = c(4, 7),  col.est = c(palette$dark, palette$red))

mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5)

Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=1)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=1)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5)
dev.off()



# Panel B - Aggregate units of analysis
fe_results <- full_summary %>% 
  filter(sigma_p==0.3,
         !((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         !pixel.fe==T) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))


coverage <- fe_results$cover
coverage[is.na(coverage)] <- 0
c_print <- f(round(fe_results$cover, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- fe_results$RMSE
RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(fe_results$RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "


schart_results <- fe_results %>% 
  select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q05, q95, pixel, county, grid, property, treatment.fe, county.fe, grid.fe, property.fe))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list(
  "Unit of analysis:" = c("pixel", "county", "grid", "property"),
  "Fixed effects:" = c("treatment FE", "grid FE", "property FE", "county FE"))


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"prop_spec.png"), width = 9, height = 7, units = "in", res = 150)
par(bg = palette$white)

schart(schart_results, labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = c(4, 7), leftmargin = 12, col.est = c(palette$dark, palette$red))

mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5)

Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=1)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=1)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5)
dev.off()



