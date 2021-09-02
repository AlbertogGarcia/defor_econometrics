library(ggplot2)
source(here::here('unbiased_dgp', 'multi_gridsize.R'))

# we start with our base parameterization without property level perturbations
std_a = 0.1
std_v = 0.5
std_p = 0.25
std_b3 = 0.25
years = 3
nobs = 10000
n = 100

ppoints = 100

# here are the landscape characteristics in this parameterization
# note that the bias for the TWFE model will be equal to the pre-treatment difference in deforestation rtes, which is 0.03
base_0 = .02
base_1 = .05
trend = -.005
ATT = -.01

seq <- seq.int(from = 1, to = 25, by = 1)
cellsize_list <- as.list(seq)

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)

cpoints = 30
set.seed(2)
# runs simulation for all models
gridsizes <- multi_gridsize(n, nobs, years, b0, b1, b2_0, b2_1, std_a, std_v, std_p, std_b3, given_ATT = ATT, cellsize_list, ppoints, cpoints, rm.selection = TRUE)
#set long dataframe with summary for each iteration
summary_grid <- gridsizes$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
            q25 = quantile(bias, probs = .25),
            q75 = quantile(bias, probs = .75),
            Bias = mean(bias),
            cover = mean(cover),
            parea = mean(parea))

grid_area_plot <- ggplot(summary_grid, aes(x = avg_garea, y = Bias))+
  geom_vline(xintercept = mean(summary_grid$parea), linetype = "dashed", color = "red") +
  #geom_errorbar(aes(ymin=q25, ymax=q75), width=.1)+
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey90")  +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed") +
  geom_line() +
  #geom_point(size=1.5, shape=21, fill="white")+ 
  ylab("Bias")+xlab("Grid area")+#+xlim(0, 625)+
  annotate(geom="text", x=301, y=0.025, label="average property area",
           color="red")+
  geom_segment(aes(x = 300, y = 0.023, xend = 105, yend = 0.015), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  theme_minimal()

ggsave(filename = "grid_size.png" , path = "paper/figs", width = 4, height = 5)

####################################################################################################

base_0 = .05
base_1 = .02
trend = .005
ATT = .01

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)


# runs simulation for all models
gridsizes_b <- multi_gridsize(n, nobs, years, b0, b1, b2_0, b2_1, std_a, std_v, std_p, std_b3, given_ATT = ATT, cellsize_list, ppoints, cpoints, rm.selection = TRUE)
#set long dataframe with summary for each iteration
summary_grid_b <- gridsizes_b$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))

grid_area_plot_b <- ggplot(summary_grid_b, aes(x = avg_garea, y = Bias))+
  geom_vline(xintercept = mean(summary_grid_b$parea), linetype = "dashed", color = "red") +
  #geom_errorbar(aes(ymin=q25, ymax=q75), width=.1)+
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey90")  +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed") +
  geom_line() +
  #geom_point(size=1.5, shape=21, fill="white")+ 
  ylab("Bias")+xlab("Grid area")+#+xlim(0, 625)+
  annotate(geom="text", x=301, y=0.025, label="average property area",
           color="red")+
  geom_segment(aes(x = 300, y = 0.023, xend = 105, yend = 0.015), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  theme_minimal()

####################################################################################################

base_0 = .04
base_1 = .03
trend = .01
ATT = -0.02

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)


# runs simulation for all models
gridsizes_c <- multi_gridsize(n, nobs, years, b0, b1, b2_0, b2_1, std_a, std_v, std_p, std_b3, given_ATT = ATT, cellsize_list, ppoints, cpoints, rm.selection = TRUE)
#set long dataframe with summary for each iteration
summary_grid_c <- gridsizes_c$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))

grid_area_plot_c <- ggplot(summary_grid_c, aes(x = avg_garea, y = Bias))+
  geom_vline(xintercept = mean(summary_grid_c$parea), linetype = "dashed", color = "red") +
  #geom_errorbar(aes(ymin=q25, ymax=q75), width=.1)+
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey90")  +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed") +
  geom_line() +
  #geom_point(size=1.5, shape=21, fill="white")+ 
  ylab("Bias")+xlab("Grid area")+#+xlim(0, 625)+
  annotate(geom="text", x=301, y=0.025, label="average property area",
           color="red")+
  geom_segment(aes(x = 300, y = 0.023, xend = 105, yend = 0.015), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  theme_minimal()

####################################################################################################

base_0 = .03
base_1 = .04
trend = -.01
ATT = 0.02

# we'll need to compute the parameters 
std_avp = (std_a^2+std_v^2+std_p^2)^.5
std_avpt = (std_a^2+std_v^2+std_p^2+std_b3^2)^.5
b0 = qnorm(base_0, mean = 0, sd = std_avp)
b1 = qnorm(base_1, mean = 0, sd = std_avp) - b0
b2_0 = qnorm(trend + base_0, mean = 0, sd = std_avp) - b0
b2_1 = qnorm(trend + base_1, mean = 0, sd = std_avp) - b0 - b1
b3_mu = qnorm( pnorm(b0+b1+b2_1, mean = 0, sd = std_avp) + ATT , mean = 0, sd = std_avpt) - (b0 + b1 + b2_1)


# runs simulation for all models
gridsizes_d <- multi_gridsize(n, nobs, years, b0, b1, b2_0, b2_1, std_a, std_v, std_p, std_b3, given_ATT = ATT, cellsize_list, ppoints, cpoints, rm.selection = TRUE)
#set long dataframe with summary for each iteration
summary_grid_d <- gridsizes_d$summary_long %>%
  mutate_at(vars(bias, cover), as.numeric)%>%
  dplyr::group_by(avg_garea)%>%
  dplyr::summarise(RMSE = rmse(bias, 0),
                   q25 = quantile(bias, probs = .25),
                   q75 = quantile(bias, probs = .75),
                   Bias = mean(bias),
                   cover = mean(cover),
                   parea = mean(parea))

grid_area_plot_d <- ggplot(summary_grid_d, aes(x = avg_garea, y = Bias))+
  geom_vline(xintercept = mean(summary_grid_d$parea), linetype = "dashed", color = "red") +
  #geom_errorbar(aes(ymin=q25, ymax=q75), width=.1)+
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "grey90")  +
  geom_hline(yintercept = 0, color = "green", linetype = "dashed") +
  geom_line() +
  #geom_point(size=1.5, shape=21, fill="white")+ 
  ylab("Bias")+xlab("Grid area")+#+xlim(0, 625)+
  annotate(geom="text", x=301, y=0.025, label="average property area",
           color="red")+
  geom_segment(aes(x = 300, y = 0.023, xend = 105, yend = 0.015), color = "red",
               arrow = arrow(length = unit(0.5, "cm")))+
  theme_minimal()
