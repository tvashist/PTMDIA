require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)

#See if there are differences in iRTs for heavies detected in both

pro_combined_lib <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_3SS_Spiked/LibraryExport/PTMDIAProject_3DIASpikedPhosphoBG_5Trans0.75LocFilter.tsv")
pro_3SS_lib <- read_tsv("")