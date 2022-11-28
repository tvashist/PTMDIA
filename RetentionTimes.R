require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)




pro_combined_lib <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_Combined/LibraryExport/Annotated_PTMDIAProject_TimsTOFPro_Combined_Library.tsv")
pro_ss_lib <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_3SS_Spiked/LibraryExport/PTMDIAProject_3DIASpikedPhosphoBG_5Trans0.75LocFilter.tsv")
# 
# heavies_report <- read_tsv("Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv")
# View(heavies_report$LibraryHeavies) 
# heavies <- na.omit(heavies_report$LibraryHeavies)
# View(heavies)
# 
# heavy_ss <- pro_ss_lib %>% filter(IntLabeledPeptide %in% heavies)
# heavy_combined <- pro_combined_lib %>% filter(IntLabeledPeptide %in% heavies)
# View(heavy_ss)
# View(heavy_combined)

combined_fractionated <- split(pro_combined_lib, pro_combined_lib$Source)[[1]]
View(combined_fractionated)

combined_SS <- split(pro_combined_lib, pro_combined_lib$Source)[[2]]
View(combined_SS)


# grouped_ss <- pro_ss_lib %>% group_by(IntLabeledPeptide) %>% summarise(PrecursorMz = PrecursorMz, iRT = iRT) %>% unique(.) 
# grouped_combined <- pro_combined_lib %>% group_by(IntLabeledPeptide) %>% summarise(PrecursorMz = PrecursorMz, iRT = iRT, ReferenceRun = ReferenceRun) %>% unique(.) %>% filter(IntLabeledPeptide %in% grouped_ss$IntLabeledPeptide)


###Comparing between IDs from 12fxn and 3SS in combined library###
grouped_ss <- combined_SS %>% group_by(IntLabeledPeptide) %>% summarise(iRT = iRT) %>% unique(.) 
grouped_fractionated <- combined_fractionated %>% group_by(IntLabeledPeptide) %>% summarise(iRT = iRT) %>% unique(.) %>% filter(IntLabeledPeptide %in% grouped_ss$IntLabeledPeptide)

sub_ss <- grouped_ss %>% filter(IntLabeledPeptide %in% grouped_fractionated$IntLabeledPeptide)



View(grouped_ss)
View(grouped_fractionated)
View(sub_ss)


rt <- ggplot(data = data.frame(x=sub_ss$iRT, y=grouped_fractionated$iRT), aes(x =x, y= y)) +
  geom_point()+
  xlab("SS_iRT") +
  ylab("Fractionated_iRT")
  
rt


mz <- ggplot(data = data.frame(x=grouped_ss$PrecursorMz, y=grouped_combined$PrecursorMz), aes(x =x, y= y)) +
  geom_point()+
  xlab("Heavies_SS_PrecMZ") +
  ylab("Heavies_Combined_PrecMZ")

mz
  

  
