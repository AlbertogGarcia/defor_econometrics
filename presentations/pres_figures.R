
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
source(here::here("unbiased_dgp", "full_landscape.R"))
library(ggpubr)

set.seed(5597)


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

results_dir <- paste0(getwd()[1], '/paper/')
out_dir <- paste0(getwd()[1], '/presentations/figs/')


base_0 = .04
base_1 = .08
trend = -.005
ATT = -.05

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
         y_counterfactual = y_counterfactual %>% recode("-1" = "Previously deforested", 
                                                        "0" = "Not deforested", 
                                                        "1" = "Deforested"),
         y = y %>% recode("-1" = "Previously deforested", 
                          "0" = "Not deforested", 
                          "1" = "Deforested"),
         treat = as.factor(treat),
         treat = treat %>% recode("0" = "Stable forest - not treated",
                                  "1" = "Stable forest - treated")) %>% 
  st_as_sf() %>%
  dplyr::select(pixels, year, treat, defor, y_it, defor_year, y, y_counterfactual) 



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Landscape plots --------------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fills = c("Previously deforested" = palette$light_grey,
          "Not deforested" = palette$dark, 
          "Deforested" = palette$red,
          "Stable forest - treated" = palette$dark,
          "Stable forest - not treated" = palette$blue)


y_fills = c("2011" = "#d0d1e6",
            "2012" = "#a6bddb",
            "2013" = "#74a9cf",
            "2014" = "#3690c0",
            "2015" = "#0570b0",
            "2016" = palette$dark,
            "Not deforested" = palette$light_grey)


# Pre-deforestation landscape
initial_df <- data_df %>%
  filter(year == 2) %>% 
  mutate(y_it = y_it %>% recode("Deforested" = "Not deforested"))

initial_forest <- ggplot() + 
  geom_sf(data = initial_df, aes(fill = y_it), color = "white", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values= fills)
initial_forest <- format_fig(initial_forest)

initial_forest %>% save_separated("initial_forest")


# Year of deforestation
plot_df <- data_df %>%
  filter(year == 2) %>% 
  mutate(defor_year = as.character(defor_year + 2010) %>%  
           recode("2017" = "Not deforested"))

defor_year <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = defor_year), color = "white", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values= y_fills)

defor_year <- format_fig(defor_year)
defor_year
defor_year %>% save_separated("defor_year")


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

# Counterfactual in year 2
i = 5
plot_df <- data_df %>%
  filter(year == i)
defor_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = y), color = "white", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values= fills)
defor_plot %>% 
  format_fig %>% 
  save_separated(paste0("only_defor_baseline", as.character(i)))
defor_plot

defor_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = y_counterfactual), color = "white", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values= fills)
defor_plot %>% 
  format_fig %>% 
  save_separated(paste0("only_defor_cf", as.character(i)))
defor_plot


# Visualizing exposure to treatment
data_df <- data_df %>% 
  mutate(all_forest = "Stable forest - treated")

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
plot_df <- data_df %>% 
  filter(year==3) %>% 
  mutate(plot_var = ifelse(y_it %in% c("Deforested", "Previously deforested"), as.character(y_it), as.character(treat)))

treat_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = intervention_area_merge, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
treat_plot %>% 
  format_fig() %>% 
  save_separated(paste0("treatment")) 


treat_plot_no_defor <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = treat), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = intervention_area_merge, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
treat_plot_no_defor %>% 
  format_fig() %>% 
  save_separated(paste0("treatment_nodefor")) 


# County-level analysis
county_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = c_bounds, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
county_plot %>% 
  format_fig() %>% 
  save_separated(paste0("county")) 

county_plot <- ggplot() + 
  geom_sf(data = initial_df, aes(fill = y_it), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = c_bounds, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
county_plot %>% 
  format_fig() %>% 
  save_separated(paste0("county_all_forest")) 


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


grid_plot <- ggplot() + 
  geom_sf(data = initial_df, aes(fill = y_it), color = "black", shape = 22, alpha = 1, size = 2.8) +
  scale_fill_manual(values = fills)  + 
  geom_vline(xintercept = c(0,8,16,24,32), color="white", size=1.25) + 
  geom_hline(yintercept = c(0,8,16,24,32), color="white", size=1.25) +
  xlab("") +
  ylab("")

grid_plot <- grid_plot %>% 
  format_fig() %>% 
  save_separated(paste0("grid_all_forest")) 


# Property-level analysis
prop_plot <- ggplot() + 
  geom_sf(data = plot_df, aes(fill = plot_var), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = p_bounds, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
prop_plot %>% 
  format_fig() %>% 
  save_separated(paste0("property")) 

prop_plot <- ggplot() + 
  geom_sf(data = initial_df, aes(fill = y_it), color = "black", shape = 22, alpha = 1, size = 2.8) +
  geom_sf(data = p_bounds, size = 1.5, fill = "NA", color = palette$white) +
  scale_fill_manual(values = fills)
prop_plot %>% 
  format_fig() %>% 
  save_separated(paste0("property_all_forest")) 



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
twfe_df <- read.csv(paste0(results_dir, "results/TWFE_long.csv")) %>% 
  select(-X)%>%
  rename(parameterization = 12)


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
  annotate("text", x = -0.002, y = 25, color = palette$white, label = "DID", size = 3) +
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
  annotate("text", x = 0.022, y = 25, color = palette$dark, label = "TWFE", size = 3) +
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
  annotate("text", x = 0.033, y = 25, color = palette$dark, label = "TWFE", size = 3) +
  annotate("text", x = 0.04, y = 210, label = "Treated baseline = 0.05\nUntreated baseline = 0.02\nATT = -0.01\nTrend = -0.005", size = 3) +
  geom_vline(xintercept = 0, color = palette$red, linetype = "longdash", size = 0.5) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  xlab("Bias") +
  ylab("Frequency")
ggsave(paste0(out_dir, "twfe_4.png"), width = 5, height = 2.6, units = "in")


plot_cox <- readRDS(paste0(results_dir, "results/summary_long.RDS"))%>%
  filter(std_p==0,
         HE.estimator == 0,
         is.na(notes),
         (pixel.fe == 0 & county.fe == 0 & property.fe == 0 & grid.fe == 0)
  ) %>%
  mutate(model = ifelse(cox == 1, "Cox PH DID", "DID"))
ggplot() +
  geom_density(data = plot_cox, aes(x=as.numeric(bias), fill = model)) +
  xlim(c(-0.008, 0.018)) +
  ylim(c(0, 250)) +
  scale_fill_manual(values = c("DID"=palette$dark, "Cox PH DID"=palette$light_grey)) +
  theme_bw(base_size = 10) +
  annotate("text", x = -0.0015, y = 25, color = palette$white, label = "DID", size = 4) +
  annotate("text", x = 0.0075, y = 30, color = palette$dark, label = "Cox PH DID", size = 4) +
  geom_vline(xintercept = 0, color = palette$red, linetype = "longdash", size = 0.5) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  xlab("Bias") +
  ylab("Frequency")  
ggsave(paste0(out_dir, "cox_bias.png"), width = 5, height = 3, units = "in")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Spec charts: Emphasize aggregation as solution -------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

full_summary <- readRDS(paste0(results_dir, "results/summary_long.RDS"))

full_summary <- full_summary %>% 
  mutate(bias = as.numeric(bias),
         cover = as.numeric(cover)) %>% 
  group_by(across(c(-iteration, -bias, -cover))) %>% 
  mutate(glabel = cur_group_id())


full_summary <- full_summary %>%
  summarize(mean_bias = mean(bias),
            cover = mean(cover),
            q025 = quantile(bias, 0.025),
            q975 = quantile(bias, 0.975),
            RMSE = sqrt(mean(bias**2)),
            glabel = cur_group_id())

## Option 1 - aggregate FEs
select_results <- full_summary %>% 
  filter(std_p==0,
         pixel == 1,
         (treatment.fe == 1 | grid.fe == 1  | property.fe == 1 | county.fe == 1),
         gridsize == 10 | is.na(gridsize),
         cox == 0,
         !(treatment.fe == 1 & property.fe == 1),
         !(property.fe == 1 & se_county == 1),
         is.na(notes)) %>% 
  mutate(cover = as.numeric(cover),
         cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel.fe), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))


# ggplot(select_results, aes(y = bias, group = glabel)) +
#   geom_density()

ggplot(select_results, aes(x = glabel, y=mean_bias)) + 
  geom_errorbar(aes(ymin=q025, ymax=q975), width=.1) +
  geom_point() +
  theme_bw()


na_row <- c("median_bias" = NA,
            "q025" = NA,
            "q975" = NA)

select_results <- rbind(na_row, select_results, na_row)%>%
  as.data.frame()


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
  select(c(mean_bias, q025, q975, pixel.fe, treatment.fe, county.fe, grid.fe, property.fe))%>%
  mutate_at(vars(pixel.fe:property.fe), ~as.logical(as.integer(.)))
index.ci <- match(c("q025","q975"), names(schart_results))

labels <- list("Pixel FE", "Treatment FE", "County FE", "Grid FE", "Property FE")


topline = -0.0075
midline = topline-0.015
ylim <- c(midline-0.012,0.0075)


png(paste0(out_dir,"spec_did_fe.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, 
       labels = labels, 
       ylim = ylim, 
       axes = FALSE, 
       index.ci=index.ci, 
       ylab="", 
       #       highlight = 2, 
       leftmargin = 15, 
       col.est = c(palette$dark, palette$red), 
       heights = c(4,2.75), 
       cex = c(1.5,1.5))

mtext("Bias   ", side=2, at = 0.00, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.005, 0, 0.005), labels=TRUE, cex = 3)
#abline(h=0, lty=2, lwd = 2)
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
  filter(std_p==0,
         pixel==0,
         weights == 0,
         gridsize == 10 | is.na(gridsize)
  ) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(pixel.fe), desc(treatment.fe), desc(county.fe), desc(grid.fe), desc(property.fe))


select_results <- rbind(na_row, select_results, na_row)%>%
  as.data.frame()


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
  select(c(mean_bias, q025, q975, county.fe, grid.fe, property.fe))%>%
  mutate_at(vars(county.fe:property.fe), ~as.logical(as.integer(.)))
index.ci <- match(c("q025","q975"), names(schart_results))

labels <- list("Grid", "County", "Property")


# topline = -0.01
# midline = topline-0.015
# ylim <- c(midline-0.012,0.01)
# 

png(paste0(out_dir,"spec_did_agg.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, 
       labels = labels, 
       ylim = ylim, 
       axes = FALSE, 
       index.ci=index.ci, 
       ylab="", 
       #      highlight = 2, 
       leftmargin = 15, 
       col.est = c(palette$dark, palette$red), 
       heights = c(4, 2.75), 
       cex = c(1.5,1.5))

mtext("Bias   ", side=2, at = 0.0, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.005, 0, 0.005), labels=TRUE, cex = 3)
#abline(h=0, lty=2, lwd = 2)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0075, paste0(RMSE_print[i]), col="black", font=1, cex=2)
  text(x= i, y=midline-0.0075, paste0(c_print[i]), col="black", font=1, cex=2)
})

mtext("RMSE", side=2, at = midline+0.0075, font=2, las=1, line=.5, cex = 2)

mtext("Coverage\nprobability", side=2, at = midline-0.0075, font=2, las=1, line=.5, cex = 2)
dev.off()



## Cox proportional hazards and proposed survival estimator
select_results <- full_summary %>% 
  filter(std_p==0,
         cox == 1
  ) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>%
  arrange(HE.estimator)

select_results <- rbind(na_row, select_results, na_row)%>%
  as.data.frame()


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
  mutate(CoxDID = ifelse(HE.estimator == 0, TRUE, FALSE),
         Alternative = ifelse(HE.estimator == 1, TRUE, FALSE))%>%
  select(c(mean_bias, q025, q975, CoxDID, Alternative))
index.ci <- match(c("q025","q975"), names(schart_results))

labels <- list("Cox PH DID", "Proposed alternative")


topline = -0.01
midline = topline-0.01
ylim <- c(midline-0.008,0.011)


png(paste0(out_dir,"cox_did.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1.5,0,1,1))

schart(schart_results, 
       labels = labels, 
       ylim = ylim, 
       axes = FALSE, 
       index.ci=index.ci, 
       ylab="", 
       highlight = 2, 
       leftmargin = 20, 
       col.est = c(palette$dark, palette$red), 
       heights = c(3, 2), 
       cex = c(1.25,1.25))

mtext("Bias", side=2, at = 0.0, font=2, las=1, line=.5, cex = 2)

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
# Spec charts: Problem of property-level disturbances ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_summary <- readRDS(paste0(results_dir, "results/summary_full.RDS"))

full_summary <- full_summary %>% 
  mutate(bias = as.numeric(bias),
         cover = as.numeric(cover)) %>% 
  filter(prop_concern == 0)%>%
  group_by(across(c(-iteration, -bias, -cover))) %>% 
  mutate(glabel = cur_group_id())


full_summary <- full_summary %>%
  summarize(mean_bias = mean(bias),
            cover = mean(cover),
            q025 = quantile(bias, 0.025),
            q975 = quantile(bias, 0.975),
            RMSE = sqrt(mean(bias**2)),
            glabel = cur_group_id())



select_results <- full_summary %>% 
  filter(treatment.fe==1,
         !((se_grid==1 | se_property==1 | se_county==1)),
         cox == 0,
         is.na(notes)) %>% 
  mutate(s0 = std_p==0,
         s1 = std_p==0.1,
         s2 = std_p==0.2,
         s3 = std_p==0.3) %>% 
  arrange(desc(std_p))
# arrange(desc(s0), desc(s1), desc(s2), desc(s3))

na_row <- c("mean_bias" = NA,
            "q025" = NA,
            "q975" = NA)

select_results <- rbind(na_row, select_results, na_row)%>%
  as.data.frame()

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
  select(c(mean_bias, q025, q975, s0, s1, s2, s3))%>%
  mutate_at(vars(s0:s3), ~as.logical(as.integer(.)))
index.ci <- match(c("q025","q975"), names(schart_results))

labels <- list("Value of sigma" = c("0", "0.1", "0.2", "0.3"))


topline = -0.015
midline = topline-0.015
ylim <- c(midline-0.012,0.006)


png(paste0(out_dir,"spec_prop_did.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, labels = labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       leftmargin = 15, col.est = c(palette$dark, palette$red), heights = c(4.25,3.5), cex = c(1.5,1.5))

mtext("Bias   ", side=2, at = 0.0, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.01, 0.005), labels=TRUE, cex = 3)
#abline(h=0, lty=2, lwd = 2)
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
  filter(std_p==0.3,
         pixel == 1,
         (treatment.fe == 1 | grid.fe == 1  | property.fe == 1 | county.fe == 1),
         gridsize == 10 | is.na(gridsize),
         cox == 0,
         !(treatment.fe == 1 & property.fe == 1),
         !(property.fe == 1 & se_county == 1),
         is.na(notes)) %>% 
  mutate(cover = as.numeric(cover),
         cover_dif = (cover - 0.95)*100) %>% 
  arrange(desc(treatment.fe), property.fe, grid.fe)

na_row <- c("mean_bias" = NA,
            "q025" = NA,
            "q975" = NA)

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
  select(c(mean_bias, q025, q975, treatment.fe, county.fe, grid.fe, property.fe))%>%
  mutate_at(vars(treatment.fe:property.fe), ~as.logical(as.integer(.)))%>%
  as.data.frame()
index.ci <- match(c("q025","q975"), names(schart_results))

labels <- list("Treatment FE", "County FE", "Grid FE", "Property FE")


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"spec_prop_fe.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)
par(oma=c(1,0,1,1))

schart(schart_results, labels = labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = 5, leftmargin = 15, col.est = c(palette$dark, palette$red), heights = c(4, 2.75), cex = c(1.5,1.5))

mtext("Bias   ", side=2, at = 0.0, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.01, 0.01), labels=TRUE, cex = 3)
#abline(h=0, lty=2, lwd = 2)
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
  filter(std_p==0.3,
         #!((se_grid==1 | se_property==1 | se_county==1) & (treatment.fe==1)),
         pixel==0,
         weights == 0,
         gridsize == 10 | is.na(gridsize)
  ) %>% 
  mutate(cover_dif = (cover - 0.95)*100) %>% 
  arrange(property.fe, grid.fe, county.fe)

na_row <- c("mean_bias" = NA,
            "q025" = NA,
            "q975" = NA)

select_results <- rbind(na_row, select_results, na_row)


coverage <- select_results$cover
# coverage[is.na(coverage)] <- 0
c_print <- f(round(coverage, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- select_results$RMSE
# RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)
RMSE_print<- f(round(RMSE, digits =4))
RMSE_print[is.na(RMSE_print)] <- " "

schart_results <- select_results %>% 
  # select(-c(se_pixel, se_grid, se_property, se_county, RMSE, cover, sigma_p, cover_dif, pixel.fe)) %>% 
  select(c(mean_bias, q025, q975, county, grid, property))%>%
  mutate_at(vars(county:property), ~as.logical(as.integer(.)))%>%
  as.data.frame()
index.ci <- match(c("q025","q975"), names(schart_results))

labels <- list("Grid", "County", "Property")


topline = -0.03
midline = topline-0.015
ylim <- c(midline-0.012,0.02)


png(paste0(out_dir,"spec_prop_agg.png"), width = 10, height = 9, units = "in", res = 150)
par(bg = palette$white)

schart(schart_results, labels, ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = 4, leftmargin = 18, col.est = c(palette$dark, palette$red), heights = c(4, 2.75), cex = c(1.5, 1.5))

mtext("Bias", side=2, at = 0.018, font=2, las=1, line=.5, cex = 2)

Axis(side=2, at = c(-0.01, 0, 0.01), labels=TRUE, cex = 3)
#abline(h=0, lty=2, lwd = 2)
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
# Multi GT results ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
county_es <- readRDS(paste0(results_dir,"results_multi/county_es.rds"))
pixel_es <- readRDS(paste0(results_dir,"results_multi/pixel_es.rds"))%>%
  filter(estimator != "Borusyak, Jaravel, Spiess (2021)")


color_scale = c("TWFE" = palette$red, "Gardner (2021)" = palette$dark, "Callaway and Sant'Anna (2020)" = palette$blue, "Borusyak, Jaravel, Spiess (2021)" = palette$dark_green,  "Truth" = palette$green)

# Get list of estimators
estimators = unique(county_es$estimator)

# Subset factor levels
levels = c("TWFE", "Borusyak, Jaravel, Spiess (2021)", "Callaway and Sant'Anna (2020)", "Gardner (2021)", "Truth")
levels = levels[levels %in% estimators]

# Subset color scales
color_scale = color_scale[names(color_scale) %in% estimators]

# create confidence intervals
out_pix = pixel_es %>%
  dplyr::mutate(
    ci_lower = q05,
    ci_upper = q95,
    estimator = factor(estimator, levels = levels)
  )

# create confidence intervals
out_county = county_es %>%
  dplyr::mutate(
    ci_lower = q05,
    ci_upper = q95,
    estimator = factor(estimator, levels = levels)
  )

# position 
position = position_dodge(width = 0.5)

ylimmin = min(min(out_pix$ci_lower), min(out_county$ci_lower))
ylimmax = max(max(out_pix$ci_upper), max(out_county$ci_upper))
ylim = c(ylimmin, ylimmax)

par(oma=c(0,0,0,0))

p1 <- ggplot(out_pix, ggplot2::aes(x = term, y = estimate, color = estimator, ymin = ci_lower, ymax = ci_upper)) +
  geom_point(position = position, size = 2.6) +
  geom_errorbar(position = position) +
  geom_vline(xintercept = -0.5, linetype = "dashed", size = 0.25) +
  labs(y = "Mean point estimate", x = "Event Time", color = "Estimator") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = color_scale) +
  guides(
    color = ggplot2::guide_legend(title.position = "top", nrow = 2)
  ) +
  ylim(ylimmin, ylimmax)+
  ggtitle(" Pixel as unit of analysis")+
  theme(plot.title = element_text(hjust = 0))
p1 + theme(legend.position = "bottom")

ggsave(paste0(out_dir, "new_did_pix.png"), width = 8, height = 6, units = "in")


p2 <- ggplot(out_county, ggplot2::aes(x = term, y = estimate, color = estimator, ymin = ci_lower, ymax = ci_upper)) +
  geom_point(position = position, size = 2.6) +
  geom_errorbar(position = position) +
  geom_vline(xintercept = -0.5, linetype = "dashed", size = 0.25) +
  labs(y = "Mean point estimate", x = "Event Time", color = "Estimator") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = color_scale) +
  guides(
    color = ggplot2::guide_legend(title.position = "top", nrow = 2)
  ) +
  ylim(ylimmin, ylimmax)+
  ggtitle("Aggregated unit of analysis")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0))
p2 + theme(legend.position = "bottom")

ggsave(paste0(out_dir, "new_did_agg.png"), width = 8, height = 6, units = "in")







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Updated figures - 9/1/2023 ------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Solution - aggregate fixed effects
specchart_long <- readRDS("paper/results/results_selection.rds")

df_summary <- specchart_long %>%
  filter(is.na(notes),
         # pixel.fe ==0,
         weights == 0 ,
         (grid.fe == 0 | gridsize == max(gridsize, na.rm = T)),
         (property.fe == 0 | se_county == 0)
  ) %>%
  mutate_at(vars(bias, cover), as.numeric
  )%>%
  mutate_at(vars(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, cox, HE.estimator, se_pixel, se_grid, se_property, se_county), ~ as.logical(as.integer(.))
  )%>%
  group_by(pixel, grid, property, county, pixel.fe, grid.fe, property.fe, county.fe, treatment.fe, cox, HE.estimator, se_pixel, se_grid, se_property, se_county)%>%
  summarise(RMSE = rmse(bias, 0),
            q05 = quantile(bias, probs = .05),
            q95 = quantile(bias, probs = .95),
            Bias = mean(bias),
            cover = mean(cover))%>%
  select(Bias, everything())

df_summary_1 <- rbind(df_summary[10,], df_summary[5,], df_summary[7:9,])
df_summary_2 <- rbind(df_summary[10,],  df_summary[1:3,])


labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "Survival" = c("Cox PH model", "ATT-Cox estimator"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))

na_row <- setNames(data.frame(matrix(ncol = ncol(df_summary), nrow = 0)), colnames(df_summary))
na_row[1,] <- NA

select_results <- rbind(na_row, as.data.frame(df_summary_1), na_row)


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
topline = -0.009

midline = topline-0.004-.0025

ylim <- c(midline-.005,0.04)
#bottomline = (min(ylim)+topline)/2
#Create the plot

rmse_cex = 1
cov_cex = 1
leg_cex = 1.2

par(oma=c(1,0,1,1))

svg(paste0('presentations/figs/1_FE_agg.svg'))
schart(select_results,labels, ylim = ylim, index.ci=index.ci, col.est = c(palette$dark, palette$red),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red),
       ylab="    Bias", highlight=c(2)
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c(0, 0.01, 0.02, 0.03), labels=TRUE)
abline(h=topline)
abline(h=midline)
#abline(h=bottomline)
lapply(1:(length(RMSE_print)), function(i) {
  # text(x= i, y=min(ylim)+.002, paste0(i), col="black", font=2, cex = 1)
  # mtext(if(i<6){paste0(i)}, side=1, at = i + 1, font=2, cex=.95)#, line=1, at=-1)
  text(x= i, y=midline+0.003, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=min(ylim)+0.0006, paste0(c_print[i]), col="black", font=1, cex=cov_cex )
})
dev.off()



## Solution - aggregate units of observation
select_results <- rbind(na_row, as.data.frame(df_summary_2), na_row)


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
topline = -0.009

midline = topline-0.004-.0025

ylim <- c(midline-.005,0.04)
#bottomline = (min(ylim)+topline)/2
#Create the plot

rmse_cex = 1
cov_cex = 1
leg_cex = 1.2

par(oma=c(1,0,1,1))

svg(paste0('presentations/figs/2_unit_agg.svg'))
schart(select_results,labels, ylim = ylim, index.ci=index.ci, col.est = c(palette$dark, palette$red),
       bg.dot=c(palette$dark, "grey95", "white", palette$red),
       col.dot=c(palette$dark, "grey95", "white", palette$red),
       ylab="    Bias", highlight=c(2)
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c(0, 0.01, 0.02, 0.03), labels=TRUE)
abline(h=topline)
abline(h=midline)
#abline(h=bottomline)
lapply(1:(length(RMSE_print)), function(i) {
  # text(x= i, y=min(ylim)+.002, paste0(i), col="black", font=2, cex = 1)
  # mtext(if(i<6){paste0(i)}, side=1, at = i + 1, font=2, cex=.95)#, line=1, at=-1)
  text(x= i, y=midline+0.003, paste0(RMSE_print[i]), col="black", font=1, cex=rmse_cex)
  text(x= i, y=min(ylim)+0.0006, paste0(c_print[i]), col="black", font=1, cex=cov_cex )
})
dev.off()


## Problem - spatial error structure
results_full <- readRDS("paper/results/results_full.rds")

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

topline = -0.0125
midline = topline-0.005
ylim <- c(midline-0.004,0.0078)

par(oma=c(1,0,1,1))

svg(paste0('presentations/figs/3_spatial_error.svg'), height = 5.5, width = 7.5)

schart(as.data.frame(schart_results), labels = labels, 
       ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="",
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
  text(x= i, y=midline+0.0025, paste0(RMSE_print[i], " "), col="black", font=1, cex=rmse_cex)
  text(x= i, y=ylim[1]+0.0015, paste0(c_print[i]), col="black", font=1, cex=cov_cex)
})
# text(x=mean(1:nrow(schart_results))
#      , y=midline-.001, "Coverage probability", col="black", font=2)
# text(x=mean(1:nrow(schart_results))
#      , y=topline-.001, "RMSE", col="black", font=2)
# mtext("RMSE", side=2, at = midline+0.003, font=2, las=1, line=.5)
# mtext("Coverage\nprobability", side=2, at = midline-0.003, font=2, las=1, line=.5)

# mtext(paste0(0), side=1, at = 2, font=2, cex=1.2, line=3)
# mtext(paste0(0.1), side=1, at = 3, font=2, cex=1.2, line=3)
# mtext(paste0(0.2), side=1, at = 4, font=2, cex=1.2, line=3)
# mtext(paste0(0.3), side=1, at = 5, font=2, cex=1.2, line=3)
# mtext(bquote("Value of " ~ sigma[p]), side=1, at = 0, font=2, cex=1.25, line=3)
dev.off()

## Solution - match scale to disturbances
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

# select_results_B <- full_summary %>% 
#   filter(pixel==0) %>% 
#   dplyr::arrange(Bias)
# # dplyr::arrange(property, county, grid, gridsize)

na_row <- setNames(data.frame(matrix(ncol = ncol(select_results_A), nrow = 0)), colnames(select_results_A))
na_row[1,] <- NA

select_results <- rbind(na_row, select_results_A, na_row)

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
  mutate(F_col = FALSE,
         F_col2 = FALSE)%>%
  select(c(Bias, F_col, F_col2, q05, q95))
index.ci <- match(c("q05","q95"), names(schart_results))

labels <- list("",
               "")

topline = -0.0125
midline = topline-0.005
ylim <- c(midline-0.004,0.0078)


par(oma=c(1,0,1,1))

svg(paste0('presentations/figs/4_spatial_agg.svg'), height = 5.5, width = 7.5)

schart(as.data.frame(schart_results), labels = labels, 
       ylim = ylim, axes = FALSE, index.ci=index.ci, ylab="", 
       highlight = c(6),  
       #  leftmargin = 16, 
       col.est = c(palette$dark, palette$red), 
       heights = c(6,1), cex = c(1.25, 1.25),
       bg.dot=c(palette$dark, "white", "white", palette$red),
       col.dot=c(palette$dark, "white", "white", palette$red)
)
Axis(side=2, at = c(-0.01, -0.005, 0, 0.005), labels=TRUE)
abline(h=topline)
abline(h=midline)
lapply(1:length(RMSE), function(i) {
  text(x= i, y=midline+0.0025, paste0(RMSE_print[i], " "), col="black", font=1, cex=rmse_cex)
  text(x= i, y=ylim[1]+0.0015, paste0(c_print[i]), col="black", font=1, cex=cov_cex)
})
dev.off()

