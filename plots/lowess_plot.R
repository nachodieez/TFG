library(ggplot2)
ggplot(data, aes(x, y)) + 
  geom_point(size=1, col = "lightskyblue3") +
  geom_smooth(method = "loess", se = FALSE, col = "green") + 
  ggtitle("Curva LOWESS") + 
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
