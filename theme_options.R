library(ggprism)

theme_options <- theme_prism(base_size = 16) + theme(
    axis.title = element_text(size = 21),
    plot.title = element_text(size = 23),
    plot.caption = element_text(size = 14, face = "italic", hjust = 0),
    legend.text = element_text(size = 18),
    legend.position = "none",
    strip.text = element_text(size = 18)
)
