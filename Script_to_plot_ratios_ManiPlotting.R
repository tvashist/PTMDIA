require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library("RColorBrewer")

setwd("Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro_SmallLibSearch_LocFilter/Pro_SmallLibSearch_LocFilter_Heavies_outputs_Found/PlottingReplicates/")


list_files = list.files(pattern = "replicates_output",recursive = T)

LL <- lapply(list_files, FUN = function(L){
  df = fread(L)
})

LL2 <- bind_rows(LL)
View(LL2)


set_conc <- 1.0



ratios_for_plotting <- LL2 %>%
  select(Spike, `Expected Ratio`) %>%
  distinct() %>%
  arrange(Spike)




ratios_for_plotting$Spike  <- factor(ratios_for_plotting$Spike, levels = ratios_for_plotting$Spike)
LL2$Spike <- factor(LL2$Spike, levels = ratios_for_plotting$Spike)


ggplot(LL2 , mapping = aes(y = 1/`Actual Ratio`, x= Spike))+
  geom_boxplot()+
  geom_jitter(aes(color=Peptide), cex = 0.5)+
  scale_y_log10()+
  theme_bw()+
  geom_point(data = ratios_for_plotting, aes(x = Spike, y = 1/`Expected Ratio`),
             shape = 95, size =10,
             color= "black", alpha = 1)+
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position = "none")


LL2 %>%
  group_by(Spike) %>%
  tally()
