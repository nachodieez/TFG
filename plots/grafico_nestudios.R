library(waffle)
parts <- data.frame(
  names = c("Affymetrix: 16", "Illumina: 7", "Agilent: 6", "RNA-seq: 8"),
  vals = c(16, 7, 6, 10)
)

waffle(parts, rows = 4, title = "")

freq = c(16, 7, 6, 8)
parts <- data.frame(
  x = rep(c("Affymetrix", "Illumina", "Agilent", "RNA-seq"), freq)
)

set.seed(123)
data <- data.frame(x = sample(LETTERS[1:6], 300, replace = TRUE))
head(data)
dim(data)

ggplot(parts, aes(x = factor(x), fill = factor(x))) +
  geom_bar() +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "black") + theme_bw() + 
  xlab("Tecnolog?a") + ylab("n? estudios")  + 
  ggtitle("N?mero de estudios inclu?dos en el trabajo divididos por tecnolog?a") +
  #scale_fill_manual(values=group.colors, name = "Tecnolog?a")
  scale_fill_discrete(name = "Tecnolog?a")

group.colors <- c(Affymetrix = "#fdf2cf", Agilent = "#dce8fb", Illumina = "#dfd6e6", `RNA-seq` = "#d8e7d5")
