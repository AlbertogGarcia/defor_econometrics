cite_df <- tibble("year" = c(2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021), 
                  "ncites" = c(2, 30, 66, 81, 111, 129, 145, 144, 193))

cite_plot <- cite_df %>% 
  ggplot(aes(x = year, y = ncites)) +
  geom_bar(stat = "identity", fill = "#0c2230") +
  ylab("Number of papers") +
  xlab("Year") +
  theme_bw(base_size = 18) +
  theme(plot.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA"),
        panel.background = element_rect(fill = "#FAFAFA",colour = "#FAFAFA")) +
  scale_x_continuous(breaks= pretty_breaks())


cite_plot
