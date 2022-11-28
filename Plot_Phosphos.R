require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)

View(All_PhosphoNumbers)

ggplot(All_PhosphoNumbers , mapping = aes(x= Library, y= `Mean Phospho Count`, fill = `Search Engine`)) +
  geom_bar(stat = "identity", position = 'dodge')+
  facet_grid(. ~ Instrument)+
  theme(axis.text.x = element_text(angle = 90))


