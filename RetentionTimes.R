require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)




pro_combined_lib <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_Combined/LibraryExport/PTMDIAProject_TimsTOFPro_Combined_Library.tsv")
pro_ss_lib <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_3SS_Spiked/LibraryExport/PTMDIAProject_3DIASpikedPhosphoBG_5Trans0.75LocFilter.tsv")

heavies <- read_tsv("Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv") %>% select (., LibraryHeavies) %>% na.omit(.)
View(heavies)



heavy_ss <- pro_ss_lib %>% filter(., IntLabeledPeptide %in% heavies)
View(heavy_ss)



