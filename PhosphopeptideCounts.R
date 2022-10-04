require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)
library(UpSetR)
library(VennDiagram)



setwd("Z:/Helium_Tan/Analyses/PhosphopeptideCounts/")

list_files = list.files(pattern = "PhosphoStatsSummary",recursive = T)

LL <- lapply(list_files, FUN = function(L){
  df = fread(L)
})



LL2 <- bind_rows(LL)
View(LL2)


p <-ggplot(data=LL2, aes(x=factor(V1+1), y=Phosphopeptides, fill = Run)) +
  geom_bar(stat = 'identity', position = 'dodge')+
  xlab('Replicate')


p









