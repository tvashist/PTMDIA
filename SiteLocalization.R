require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)

lights<- read.table('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Lights.tsv', sep ='\t', header = TRUE, quote = '')
lights <- lights$Modified


setwd('Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/')

conditions <- c(0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0)

spectronaut <- read.table('directDIA/Pro/20220714_100348_PTMDIAProject_Pro_PhosphoBG_Report.tsv', sep = '\t', header= TRUE, quote = '') %>% mutate(conc = as.numeric(gsub(R.Condition, pattern = "fmol", replacement = "")))



find_localization <- function(level){
  
  spike <- spectronaut %>% filter(conc == level) %>% filter(EG.IntPIMID %in% lights) #Dataframe with all the lights detected at that specific concentration

  ggplot(spike,aes(x=EG.PTMAssayProbability))+
    geom_histogram(binwidth = 0.05)+
    xlab("PTM Localization Probability")+
    geom_vline(aes(xintercept = mean(spike$EG.PTMAssayProbability)),col='purple',size=2)+
    geom_textvline(label = format(round(mean(spike$EG.PTMAssayProbability),4), nsmall = 2), xintercept=mean(spike$EG.PTMAssayProbability), vjust = -0.015, size = 5 )+
    ggtitle(paste(level, "fmol", sep = " "))
  
  
  
}

lapply(conditions[1:2], find_localization)





  

