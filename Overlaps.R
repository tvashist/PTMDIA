require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)
library(UpSetR)
library(ggvenn)
library(readxl)
library(VennDiagram)
install.packages("pacman")


twelvefxn_search1 <- read_tsv('PTMDIAProject_SpectralLibraries/Exports/Pro_12fxnOnly/PTMDIAProject_TimsTOFPro_12fxnOnly.xls')


setwd('Z:/Helium_Tan/')

#All DIA data files
twelvefxn_search <- read.table('PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro_12fxnOnlySearch_LocFilter/20220906_093527_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_12fxnLib0.75Loc_Report.tsv', sep = '\t', header = TRUE, quote = '')
singleshots_search <- read.table('PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro_SmallLibSearch_LocFilter/20220902_105134_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_SmallLib0.75Loc_Report.tsv', sep = '\t', header = TRUE, quote = '')
combined_search <- read.table('PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro/20220803_134456_PTMDIAProject_DIACurveAnalysis_WithSpecLib_Report_addedFGLabel.tsv', sep = '\t', header = TRUE, quote = '')

#Library export files
twelvefxn_lib <- read.table('PTMDIAProject_SpectralLibraries/Exports/Pro_12fxnOnly/PTMDIAProject_TimsTOFPro_12fxnOnly.xls', sep = '\t', header = TRUE, quote = '')
singleshots_lib <- read.table('PTMDIAProject_SpectralLibraries/Exports/Pro_3DDASpikedOnly/PTMDIAProject_3DIASpikedPhosphoBG_5Trans0.75LocFilter.xls', sep = '\t', header = TRUE, quote = '')
combined_lib <- read.table('PTMDIAProject_SpectralLibraries/Exports/Pro/PTMDIAProject_Pro_PhosphoSL.xls', sep = '\t', header = TRUE, quote = '')





make_upset <- function(DIA, Lib) {
  
  search <- unique(DIA$PEP.StrippedSequence)
  library <- unique(Lib$StrippedPeptide)
  
  
  overlap <- list(DIA = search,Lib = library)
  
  u <- upset(fromList(overlap),sets.bar.color=c("maroon","blue"))
  u
  
  perc_found <- (length(search)/length(library))*100
  print(perc_found)
  
}

make_upset(twelvefxn_search,twelvefxn_lib)





          
                             
                             
                             
                             
          
