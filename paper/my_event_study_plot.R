
my_event_study_plot = function(out, seperate = TRUE, horizon = NULL) {
  
  # Get list of estimators
  estimators = unique(out$estimator)
  
  # Subset factor levels
  levels = c("TWFE", "Borusyak, Jaravel, Spiess (2021)", "Callaway and Sant'Anna (2020)", "Gardner (2021)", "Roth and Sant'Anna (2021)",  "Sun and Abraham (2020)", "Truth")
  levels = levels[levels %in% estimators]
  
  # Subset color scales
  color_scale = c("TWFE" = "#374E55", "Gardner (2021)" = "#DF8F44", "Callaway and Sant'Anna (2020)" = "#00A1D5", "Sun and Abraham (2020)" = "#B24745", "Roth and Sant'Anna (2021)" = "#79AF97", "Borusyak, Jaravel, Spiess (2021)" = "#6A6599", "Truth" = "limegreen")
  color_scale = color_scale[names(color_scale) %in% estimators]
  
  # create confidence intervals
  out = out %>%
    dplyr::mutate(
      ci_lower = q05,
      ci_upper = q95,
      estimator = factor(estimator, levels = levels)
    )
  
  # position depending on sepreate
  if(seperate) position = "identity" else position = position_dodge(width = 0.6)
  
  # Subset plot if horizon is specified
  if(!is.null(horizon)) {
    out = out %>%
      dplyr::filter(term >= horizon[1] & term <= horizon[2])
  }
  
  # max and min of limits
  y_lims = c(min(out$ci_lower), max(out$ci_upper)) * 1.05
  x_lims = c(min(out$term) - 1, max(out$term) + 1)
  
  ggplot2::ggplot(out, ggplot2::aes(x = term, y = estimate, color = estimator, ymin = ci_lower, ymax = ci_upper)) +
    { if(seperate) ggplot2::facet_wrap(~ estimator, scales="free") } +
    ggplot2::geom_point(position = position, size = 2.6) +
    ggplot2::geom_errorbar(position = position) +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed") +
    #ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(y = "Mean point estimate", x = "Event Time", color = "Estimator") +
    { if(seperate) ggplot2::scale_y_continuous(limits = y_lims) } +
    { if(seperate) ggplot2::scale_x_continuous(limits = x_lims) } +
    ggplot2::theme_minimal(base_size = 18) +
    ggplot2::scale_color_manual(values = color_scale) +
    ggplot2::guides(
      color = ggplot2::guide_legend(title.position = "top", nrow = 2)
    ) +
    ggplot2::theme(legend.position = "bottom")
  
}