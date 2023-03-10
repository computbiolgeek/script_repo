library(tidyverse)
library(PNWColors)

# load dataset
score_vs_rms <- score_file <- read_table(
  file = "score_file.scores",
  col_names = TRUE
)

# get the top model
top_model <- score_vs_rms %>% filter(
  description == "top_model_name"
)
the_rest <- score_vs_rms %>% filter(
  irmsd <= 5
)

plot_margin <- margin(
  t = 0.5, r = 0.5, b = 0.3, l = 0.3, unit = "cm"
)
funnel_plot <- ggplot(
  data = the_rest,
  mapping = aes(x = irmsd, y = ddg)
) + 
  geom_point(
    size = 2,
    alpha = 0.5
  ) + 
  geom_point(
    data = top_model,
    mapping = aes(x = irmsd, y = ddg),
    size = 4,
    color = "red",
    alpha = 0.5
  ) +
  xlab(
    label = expression("Interface r.m.s.d. (" ~ ring(A) ~ ")")
  ) +
  scale_x_continuous(
    limits = c(0, 2)
  ) + 
  scale_y_continuous(
    name = "Binding energy (REU)",
    limits = c(-440, -360),
    breaks = seq(-440, -360, 20)
  ) +
  theme(
    panel.border = element_rect(colour = "black", size = 1, fill = "transparent"),
    axis.text = element_text(size = 24, color = "black"),
    axis.title.x = element_text(
      color = "black", size = 28, margin = margin(t = 5)
    ),
    axis.title.y = element_text(
      color = "black", size = 28, margin = margin(r = 10)
    ),
    # axis.title.x = element_blank(),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks = element_line(size = 1.0),
    legend.position = "none",
    aspect.ratio = 1.0,
    plot.margin = plot_margin
  )

ggsave(
  filename = "score_vs_rmsd_funnel_plot.png",
  plot = funnel_plot,
  width = 8,
  height = 8,
  units = "in",
  device = "png",
)
