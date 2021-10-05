source(here::here('multigroup_dgp', 'multipleGT.R'))
source(here::here('multigroup_dgp', 'multipleGT_agg.R'))
source(here::here('multigroup_dgp', 'multipleGT_pix.R'))
source(here::here('multigroup_dgp', 'my_event_study_plot.R'))
library(staggered)
base_a = .07
base_b = .05
base_c = .02
trend1 = .00
trend2 = .00
trend3 = .00
ATT = -.02
dyn_ATT = 0

std_a = 0.0
std_v = 0.5
std_p = 0.0

cpoints = 70
cellsize=10 
ppoints=75

nobs = 120^2
n=100

multiGT <- multipleGT(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT, dyn_ATT = 0, std_a = 0.0, std_v = 0.25, std_p = 0.0, cellsize=10, ppoints=50, cpoints)
  
library(rio)
export(multiGT$es_long, "es_long_main.rds")

es <- multiGT$es_long %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))
county_es <- subset(es, uoa=="county"|estimator=="Truth")%>%
  filter(estimator!="Sun and Abraham (2020)")%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)
pixel_es <- subset(es, uoa=="pixel"|estimator=="Truth")%>%
  filter(estimator!="Sun and Abraham (2020)")%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)


my_event_study_plot(pixel_es, seperate = FALSE)+
  ggtitle("estimates with pixel unit of analysis")+
  geom_segment(aes(x = -2.5, y = 0, xend = -0.5, yend = 0), color = "limegreen")+
  geom_segment(aes(x = -0.5, y = -0.02, xend = 2.5, yend = -0.02), color = "limegreen")+
  ylim(-0.04, 0.05)
my_event_study_plot(county_es, seperate = FALSE)+
  ggtitle("estimates with aggregated unit of analysis (county)")+
  geom_segment(aes(x = -2.5, y = 0, xend = -0.5, yend = 0), color = "limegreen")+
  geom_segment(aes(x = -0.5, y = -0.02, xend = 2.5, yend = -0.02), color = "limegreen")+
  ylim(-0.04, 0.05)


plot_df <- panels %>%
  group_by(G, year)%>%
  summarise(defor = mean(y, na.rm=TRUE))%>%
  mutate(Group = ifelse(G==0, "never", "early group\n(treated in year 3)"),
         Group = ifelse(G==4, "late group\n(treated in year 4", Group))


ggplot(data=plot_df, aes(x=year, y=defor, colour=Group))+
  geom_line(size=1.5)+
  ylab("deforestation rate")+
  scale_y_continuous(labels = scales::percent)+
  ggtitle("Annual deforestation with multiple groups and periods")+
  theme_minimal()+
  theme(#legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm")) +
  theme(legend.title.align = 0.5)
  
##########################################################################
##########################################################################
dyn_ATT_a = 0.01
dyn_ATT_b =  -0.02
ATT_a = -0.03
ATT_b =  -0.02
set.seed(0930)
multiGT_agg <- multipleGT_agg(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize=10, ppoints=50, cpoints)

county_es <- multiGT_agg$es_long %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  filter(estimator!="Roth and Sant'Anna (2021)")%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

my_event_study_plot(county_es, seperate = FALSE)+
  ggtitle("estimates with aggregated unit of analysis (county)")+
  geom_segment(aes(x = -2.5, y = 0, xend = -0.5, yend = 0), color = "limegreen")
  #geom_segment(aes(x = -0.5, y = -0.02, xend = 2.5, yend = -0.02), color = "limegreen")


multiGT_pix <- multipleGT_pix(n, nobs, base_a, base_b, base_c, trend1, trend2, trend3, ATT_a, ATT_b, dyn_ATT_a, dyn_ATT_b, std_a, std_v, std_p, cellsize=10, ppoints=50, cpoints)

pixel_es <- multiGT_pix$es_long %>%
  group_by(term, estimator, uoa)%>%
  mutate(estimate = as.numeric(estimate))%>%
  summarise(q05 = quantile(estimate, probs = 0.05),
            q95 = quantile(estimate, probs = 0.95),
            estimate = mean(estimate))%>%
  filter(estimator!="Roth and Sant'Anna (2021)")%>%
  mutate_at(vars(term, q05,q95, estimate), as.numeric)

my_event_study_plot(pixel_es, seperate = FALSE)+
  ggtitle("estimates with pixel unit of analysis")+
  geom_segment(aes(x = -2.5, y = 0, xend = -0.5, yend = 0), color = "limegreen")

