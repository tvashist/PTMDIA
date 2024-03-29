require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)


setwd("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/DIANN/Library_3SS_Spiked/SearchOutputs/TimsTOF_Pro_Heavies_outputs_Found")


list_files = list.files(pattern = "_output",recursive = T)

LL <- lapply(list_files, FUN = function(L){
  df = fread(L)
})




LL2 <- bind_rows(LL) %>%
  filter(!is.na(`Mean Quant`)) %>%
  mutate(conc = as.numeric(gsub(Spike, pattern = "fmol", replacement = "")))

ten <- LL2 %>% filter(Spike == '10.0fmol')
View(ten)


ratios_for_plotting <- LL2 %>%
  select(Spike, `Expected Ratio`, conc) %>%
  distinct() %>%
  arrange(conc)

LL2 %>%
  group_by(Spike) %>%
  tally()


ratios_for_plotting$Spike  <- factor(ratios_for_plotting$Spike, levels = ratios_for_plotting$Spike)
LL2$Spike <- factor(LL2$Spike, levels = ratios_for_plotting$Spike)
# View(LL2)

ggplot(LL2 , aes(y = 1/`Actual Ratios`, x= Spike, color=Spike ))+
  geom_boxplot()+
  geom_jitter(size = 0.5)+
  scale_y_log10()+
  theme_bw()+
  geom_point(data = ratios_for_plotting, aes(x = Spike, y = 1/`Expected Ratio`),
             shape = 95, size =10,
  color= "black", alpha=0.5)+
  theme(axis.text.x = element_text(angle=90))




actual_median <- c()
errors <- c()

for (point in ratios_for_plotting$conc) {
  
  sub <- LL2 %>% filter (conc == point)    #Filter down to just the rows for that spike level
  
  actual <- median(sub$`Actual Ratios`)
  actual_median <- append(actual_median, actual)
  
  expected <- unique(sub$`Expected Ratio`)
  
  #Percent error between expected and actual ratios
  
  numerator <- abs(expected - actual)
  denominator <- abs(expected)
  perc_error <- (numerator/denominator) * 100
  
  errors <- append(errors, perc_error)
  
  
}

ratios_for_plotting$`Percent Error` <- errors
ratios_for_plotting$`Actual Ratio Median` <- actual_median

ratios_for_plotting <- select(ratios_for_plotting, Spike, `Expected Ratio`, `Actual Ratio Median`, `Percent Error`)


write_tsv(ratios_for_plotting, "Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro_SmallLibSearch/Pro_SmallLibSearch_Heavies_outputs_Found/PercentError.tsv")

ggplot(data = ratios_for_plotting, aes(x = Spike, y=`Percent Error`))+
  geom_point()+
  ggtitle("Percent Error Between Expected and Actual Ratios")





