require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)
library(tidyverse)



setwd('Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_SCP/Spectronaut/Library_3SS_Spiked/SearchOutputs/TimsTOF_SCP_UnfilteredSiteLocalization/')


dist <- read_tsv('Lights_UnfilteredSiteLocalization.tsv')  %>% select(2:10) %>% gather()

perc_loc <- read_tsv('Heavy_PercLocalized.tsv') %>% select(2:10) %>% gather()





SCP_levels <- c('0.004 fmol', '0.01 fmol', '0.02 fmol', '0.04 fmol','0.1 fmol', '0.2 fmol', '0.4 fmol', '1.0 fmol', '2.0 fmol' )

other_levels <- c('0.1 fmol', '0.2 fmol', '0.4 fmol', '1.0 fmol', '2.0 fmol', '4.0 fmol', '10.0 fmol')


ggplot(dist, aes(value)) +
  geom_histogram(binwidth = 0.1) +
  facet_grid(~factor(key, levels = SCP_levels )) +
  xlab('pSTY Site Localization Score')+
  ggtitle('Spectronaut')

ggplot(perc_loc, aes(x= key, y = value)) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(limits = other_levels) +
  xlab('Spike') +
  ylab('Percent Localized Sites')+
  geom_text(aes(label = round(value,2), vjust = -0.2))



# hist(diann$PTM.Site.Confidence)



