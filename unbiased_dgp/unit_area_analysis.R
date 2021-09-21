library(ggplot2)
source(here::here('unbiased_dgp', 'multiple_gridsize.R'))
source(here::here('unbiased_dgp', 'multiple_gridsize_clean.R'))

# we start with our base parameterization without property level perturbations
years = 2
nobs = 10000
n = 250

ppoints = 100

# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

seq <- seq.int(from = 1, to = 20, by = 1)
cellsize_list <- as.list(seq)

std_a = 0
std_v = 0.5
std_p = 0


# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)
b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)

# runs simulation for all models
set.seed(2)
gridsizes <- multiple_gridsize_clean(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list, ppoints, cpoints, rm.selection = FALSE)

# create dataframe
summary_grid <- gridsizes$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea, grid, pixel)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))%>%
  mutate(parameterization = "no property influence")

####################################################################################################

std_p = 0.5

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3 = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avp) - (b0 + b1 + b2_1)
b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)


# runs simulation for all models
gridsizes_pu <- multiple_gridsize_clean(n, nobs, years, b0, b1, b2_0, b2_1, b3, std_a, std_v, std_p, cellsize_list, ppoints, cpoints, rm.selection = TRUE)
#set long dataframe with summary for each iteration
summary_grid_pu <- gridsizes_pu$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea, grid, pixel)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))%>%
  mutate(parameterization = "property unobservables")

prop_sizes <- gridsizes_pu$pixloc_df %>%
  group_by(property) %>%
  summarise(psize = mean(parea))
####################################################################################################


####################################################################################################
####################################################################################################
####################################################################################
####################################################################################################

all_gridsizes <- rbind(summary_grid, summary_grid_pu)

ua_grid <- all_gridsizes %>%
  filter(grid == 1)%>%
  mutate(`unit of analysis` = "grid")

ua_pixel <- all_gridsizes %>%
  filter(grid == 0)%>%
  mutate(`unit of analysis` = "pixel")

# all_gridsizes <- rbind(ua_pixel, ua_grid)


# all_grid_plot <- ggplot(ua_pixel, aes(x = avg_garea, y = Bias, colour = parameterization#, linetype = `unit of analysis`
#                                       ))+
#   #geom_density(data = prop_sizes)+
#   geom_vline(xintercept=c(4, 49,100, 144, 225, 400), colour="grey90") +
#   #geom_vline(xintercept = 250, linetype = "dashed", color = "red") +
#   #geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey90")  +
#   geom_hline(yintercept = 0, color = "green", linetype = "dashed", size = 1.25) +
#   geom_line(size = 1.25) +
#   #geom_point(size=1.5, shape=21, fill="white")+ 
#   ylab("Bias")+#xlab("Grid area")+
#   geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
#   ylim(-0.001, 0.006)+
#   annotate(geom="text", x=201, y=0.0035, label="average property area",
#            color="red", size = 5)+
#   geom_segment(aes(x = 200, y = 0.0033, xend = 105, yend = 0.0025), color = "red",
#                arrow = arrow(length = unit(0.5, "cm")))+
#   #+xlim(1.1, 500)+
#   scale_x_continuous(
#     breaks = c(4, 49,100, 144, 225, 400),
#     labels = c("2x2", "7x7", "10x10", "12x12", "15x15", "20x20"),
#     #limits = c(0, 20^2),
#     name = "Grid size") +
#   theme_minimal()+
#   theme(text = element_text(size=15),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.x = element_blank()) 

all_grid_plot <- ggplot()+
  geom_density(data = prop_sizes, aes(psize), color = "white", fill = "grey90")+
  geom_line(data = ua_pixel, aes(x = avg_garea, y = Bias, colour = parameterization), size = 1.25) +
  geom_vline(xintercept=c(4, 49,100, 144, 225, 400), colour="grey90") +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed", size = 1.25) +
  ylab("Bias")+#xlab("Grid area")+
  geom_vline(xintercept = mean(all_gridsizes$parea), linetype = "dashed", color = "red")+
  #ylim(-0.001, 0.01)+
  annotate(geom="text", x=201, y=0.0035, label="average property area",
           color="red", size = 5)+
  geom_segment(aes(x = 200, y = 0.0033, xend = 105, yend = 0.0025), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  #+xlim(1.1, 500)+
  scale_x_continuous(
    breaks = c(4, 49,100, 144, 225, 400),
    labels = c("2x2", "7x7", "10x10", "12x12", "15x15", "20x20"),
    #limits = c(0, 20^2),
    name = "Grid size") +
  theme_minimal()+
  theme(text = element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) 

