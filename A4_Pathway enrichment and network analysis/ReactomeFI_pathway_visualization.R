# bubble plot (pathway * size)

getwd()
setwd("D:/Analysis/Biomarker/")

Sys.glob("*.csv")

df <- as.data.frame(read.csv("ReactomeFI_marker_PEA.csv"))

library(ggplot2)
library(dplyr)

df <- df %>% mutate(logFDR = -log10(FDR))
# add to -log10 FDR

df_top <- df %>% slice_head(n = 8)

p <- ggplot(df_top, aes(x= logFDR, y= reorder(Pathway, logFDR))) +
  geom_point(aes(size = Size, color = logFDR)) +
  scale_size_continuous(range = c(3, 10)) +             # bubble size (range)
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "-log10(FDR)", y = "Pathway",
       size = "Marker", color = "Size")+
  theme(axis.text.y = element_text(size = 10))

p

ggsave("top8_pathway_plot.png", plot = p,
       width = 9.64, height = 6.3, dpi = 300, units = "in")
