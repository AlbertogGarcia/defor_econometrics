library(ggplot2)
source(here::here('unbiased_dgp', 'multiple_gridsize_clean.R'))

# we start with our base parameterization without property level perturbations
years = 3
nobs = 150^2
n = 250

ppoints = 75
cpoints = 30
# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

seq <- c(1,2, 3, seq.int(from = 4, to = 22, by = 2), seq.int(from = 24, to = 100, by = 4))
seq
cellsize_list <- as.list(seq)

std_a = 0.1
std_v = 0.5

####################################################################################################


std_p = 0
# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
#std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)

# runs simulation for all models
set.seed(0930)

gridsizes <- multiple_gridsize_clean(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list, ppoints, cpoints, rm.selection = FALSE)

# create dataframe
summary_grid <- gridsizes$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea, grid, pixel, property)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))%>%
  mutate(parameterization = "no property influence")


####################################################################################################
####################################################################################################
std_p = 0.1

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
#std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)
#b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)


# runs simulation for all models
gridsizes_1 <- multiple_gridsize_clean(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list, ppoints, cpoints, rm.selection = FALSE)
#set long dataframe with summary for each iteration
summary_grid_1 <- gridsizes_1$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea, grid, pixel, property)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))%>%
  mutate(parameterization = "property unobservables (0.1)")


prop_sizes <- gridsizes_1$pixloc_df %>%
  group_by(property) %>%
  summarise(psize = mean(parea))
library(rio)
export(prop_sizes, "prop_sizes.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

std_p = 0.25
# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
#std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)



gridsizes_25 <- multiple_gridsize_clean(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list, ppoints, cpoints, rm.selection = FALSE)
#set long dataframe with summary for each iteration
summary_grid_25 <- gridsizes_25$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea, grid, pixel, property)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))%>%
  mutate(parameterization = "property unobservables (0.25)")


####################################################################################
####################################################################################################

all_gridsizes <- rbind(summary_grid, summary_grid_1, summary_grid_25)#, summary_grid_3)
export(all_gridsizes, "all_gridsizes_70a01_3.rds") #years = 3, std_a = 0.1
#export(all_gridsizes, "all_gridsizes_300a01.rds") #years = 2, std_a = 0.1
####################################################################################################

####################################################################################################
####################################################################################################
####################################################################################
####################################################################################################

library(dplyr)
all_gridsizes <- readRDS("unbiased_dgp/results/all_gridsizes_300a01_3.rds")%>%
  mutate_at(vars(avg_garea, Bias), as.numeric)
prop_sizes <- readRDS("unbiased_dgp/results/prop_sizes.rds")

ua_grid <- all_gridsizes %>%
  filter(grid == 1)%>%
  mutate(`unit of analysis` = "grid")

ua_prop0.25 <- subset(all_gridsizes, property == "unit of analysis" & parameterization == "property unobservables (0.25)")
ua_prop0.1 <- subset(all_gridsizes, property == "unit of analysis" & parameterization == "property unobservables (0.1)")
ua_prop0.0 <- subset(all_gridsizes, property == "unit of analysis" & parameterization == "no property influence")

ua_pixel <- all_gridsizes %>%
  filter(grid == 0 & pixel == 1)%>%
  mutate(`unit of analysis` = "pixel")

fe_prop0.25 <- subset(all_gridsizes, property == "fixed effects" & parameterization == "property unobservables (0.25)")
fe_prop0.1 <- subset(all_gridsizes, property == "fixed effects" & parameterization == "property unobservables (0.1)")
fe_prop0.0 <- subset(all_gridsizes, property == "fixed effects" & parameterization == "no property influence")


library(ggplot2)

my_breaks = c(9, 25, 30*30, 50*50, 75*75, 100^2)
my_labels = c("9x9", "35x5", "30x30", "50x50", "75x75", "")
xlimits = c(9, max(ua_grid$avg_garea))

ua_grid_bias <- ggplot()+
  geom_density(data = prop_sizes, aes(psize),  fill = "grey80")+
  geom_vline(xintercept=my_breaks, colour="grey90") +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed", size = 1.25) +
  geom_hline(yintercept = c(ua_prop0.25$Bias, ua_prop0.1$Bias, ua_prop0.0$Bias), color = c("blue", "green", "red"), linetype = "dashed", size = 1) +
  geom_line(data = ua_grid, aes(x = avg_garea, y = Bias, colour = parameterization), size = 1.25) +
  ylab("Bias")+#xlab("Grid area")+
  geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
  ylim(-0.004, 0.0075)+
  annotate(geom="text", x=75*75, y=0.005, label="average property area",
           color="red", size = 5)+
  geom_segment(aes(x = 75*75, y = 0.0045, xend = mean(all_gridsizes$parea) + 10, yend = 0.004), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  scale_x_continuous(
    breaks = my_breaks,
    labels = my_labels,
    limits = xlimits,
    name = "Grid size") +
  theme_minimal()+
  theme(text = element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  ggtitle("bias vs. grid size (as unit of analysis)")

ggsave(filename = "gridsize_bias_uoa.png", height = 6, width = 12)


ua_pixel_bias <- ggplot()+
  #geom_ribbon(data = ua_pixel, aes(x=avg_garea, ymin=q25, ymax=q75, fill = parameterization), alpha = 0.1)+
  geom_density(data = prop_sizes, aes(psize), fill = "grey80")+
  geom_vline(xintercept=my_breaks, colour="grey90") +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed", size = 1.25) +
  #geom_hline(yintercept = c(fe_prop0.25$Bias, fe_prop0.1$Bias, fe_prop0.0$Bias), color = c("blue", "green", "red"), linetype = "dashed", size = 1) +
  geom_line(data = ua_pixel, aes(x = avg_garea, y = Bias, colour = parameterization), size = 1.25) +
  ylab("Bias")+#xlab("Grid area")+
  geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
  ylim(-0.004, 0.0075)+
  annotate(geom="text", x=75*75, y=0.005, label="average property area",
           color="red", size = 5)+
  geom_segment(aes(x = 75*75, y = 0.0045, xend = mean(all_gridsizes$parea) + 10, yend = 0.004), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  #+xlim(1.1, 500)+
  scale_x_continuous(
    breaks = my_breaks,
    labels = my_labels,
    limits = xlimits,
    name = "Grid size") +
  theme_minimal()+
  theme(text = element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  ggtitle("bias vs. grid size (as fixed effects)")

ggsave(filename = "gridsize_bias_fe.png", height = 6, width = 12)


ua_grid_cover <- ggplot()+
  geom_density(data = prop_sizes, aes(psize), fill = "grey80")+
  geom_vline(xintercept=my_breaks, colour="grey90") +
  geom_hline(yintercept = 0.95, color = "green", linetype = "dashed", size = 1) +
  geom_line(data = ua_grid, aes(x = avg_garea, y = cover, colour = parameterization), size = 1.25) +
  ylab("Bias")+#xlab("Grid area")+
  geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
  ylim(0.5, 1)+
  annotate(geom="text", x=50*50, y=0.78, label="average property area",
           color="red", size = 5)+
  geom_segment(aes(x = 50*50, y = 0.73, xend = mean(all_gridsizes$parea) + 10, yend = 0.6), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  scale_x_continuous(
    breaks = my_breaks,
    labels = my_labels,
    limits = xlimits,
    name = "Grid size") +
  theme_minimal()+
  theme(text = element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  ggtitle("coverage vs. grid size (as unit of analysis)")

ggsave(filename = "gridsize_cover_uoa.png", height = 6, width = 12)


ua_pixel_cover <- ggplot()+
  geom_density(data = prop_sizes, aes(psize), fill = "grey80")+
  geom_vline(xintercept=my_breaks, colour="grey90") +
  geom_hline(yintercept = 0.95, color = "green", linetype = "dashed", size = 1) +
  geom_line(data = ua_pixel, aes(x = avg_garea, y = cover, colour = parameterization), size = 1.25) +
  ylab("Bias")+#xlab("Grid area")+
  geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
  ylim(0.5, 1)+
  annotate(geom="text", x=50*50, y=0.78, label="average property area",
           color="red", size = 5)+
  geom_segment(aes(x = 50*50, y = 0.73, xend = mean(all_gridsizes$parea) + 10, yend = 0.6), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  scale_x_continuous(
    breaks = my_breaks,
    labels = my_labels,
    limits = xlimits,
    name = "Grid size") +
  theme_minimal()+
  theme(text = element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  ggtitle("coverage vs. grid size (as fixed effects)")

ggsave(filename = "gridsize_cover_fe.png", height = 6, width = 12)


ua_grid_rmse <- ggplot()+
  geom_density(data = prop_sizes, aes(psize), fill = "grey80")+
  geom_vline(xintercept=my_breaks, colour="grey90") +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed", size = 1.25) +
  geom_line(data = ua_grid, aes(x = avg_garea, y = RMSE, colour = parameterization), size = 1.25) +
  ylab("Bias")+#xlab("Grid area")+
  geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
 # ylim(0, 0.0075)+
  annotate(geom="text", x=50*50, y=0.0018, label="average property area",
           color="red", size = 5)+
  geom_segment(aes(x = 50*50, y = 0.0013, xend = mean(all_gridsizes$parea) + 10, yend = 0.001), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  scale_x_continuous(
    breaks = my_breaks,
    labels = my_labels,
    limits = xlimits,
    name = "Grid size") +
  theme_minimal()+
  theme(text = element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  ggtitle("RMSE vs. grid size (as unit of analysis)")

ggsave(filename = "gridsize_rmse_uoa.png", height = 6, width = 12)


density <- c(
   "distribution of property areas" = "grey80"
)

ua_pixel_rmse <- ggplot()+
  geom_density(data = prop_sizes, aes(psize), fill = "grey80")+
  geom_vline(xintercept=my_breaks, colour="grey90") +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed", size = 1.25) +
  geom_line(data = ua_pixel, aes(x = avg_garea, y = RMSE, colour = parameterization), size = 1.25) +
  ylab("Bias")+#xlab("Grid area")+
  geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
  ylim(0, .0075)+
  annotate(geom="text", x=50*50, y=0.0056, label="average property area",
           color="red", size = 5)+
  geom_segment(aes(x = 50*50, y = 0.0053, xend = mean(all_gridsizes$parea) + 10, yend = 0.0045), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  #+xlim(1.1, 500)+
  scale_x_continuous(
    breaks = my_breaks,
    labels = my_labels,
    limits = xlimits,
    name = "Grid size") +
  theme_minimal()+
  theme(text = element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  ggtitle("RMSE vs. grid size (as fixed effects)")

ggsave(filename = "gridsize_rmse_fe.png", height = 6, width = 12)
