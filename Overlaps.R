require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)
library(UpSetR)
library(ggvenn)



setwd("Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/")


#DIA data
directDIA <- read.table(file = 'directDIA/Pro/20220714_100348_PTMDIAProject_Pro_PhosphoBG_Report.tsv', sep = '\t', header = TRUE)
combined <- read.table(file = 'SpectralLibSearch/Pro/20220803_134456_PTMDIAProject_DIACurveAnalysis_WithSpecLib_Report_addedFGLabel.tsv', sep = '\t', header = TRUE)
ThreeSS <- read.table(file = 'SpectralLibSearch/Pro_SmallLibSearch_LocFilter/20220902_105134_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_SmallLib0.75Loc_Report.tsv', sep = '\t', header = TRUE)
twelvefxn <- read.table(file = 'SpectralLibSearch/Pro_12fxnOnlySearch_LocFilter/20220906_093527_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_12fxnLib0.75Loc_Report.tsv', sep = '\t', head = TRUE)


dDIA_peptides <- as.list(unique(directDIA$PEP.StrippedSequence))
combined_peptides <- as.list(unique(combined$PEP.StrippedSequence))
threeSS_peptides <- as.list(unique(ThreeSS$PEP.StrippedSequence))
twelvefxn_peptides <- as.list(unique(twelvefxn$PEP.StrippedSequence))



DIApeptideSets <- list(
  directDIA = dDIA_peptides,
  Combined = combined_peptides,
  SingleShot = threeSS_peptides,
  TwelveFxn = twelvefxn_peptides
)

upset(fromList(peptideSets), order.by = "freq")


#Read in spectral libraries
setwd("Z:/Helium_Tan/PTMDIAProject_SpectralLibraries/Exports/")

Lib_combined <- read.table(file = 'Pro/PTMDIAProject_Pro_PhosphoSL.tsv', sep = '\t', header = TRUE)
Lib_12fxn <- read.table(file = 'Pro_12fxnOnly/PTMDIAProject_TimsTOFPro_12fxnOnly.tsv', sep = '\t', header = TRUE)

Lib_3SS <- read.table(file = 'Pro_3DDASpikedOnly/PTMDIAProject_3DIASpikedPhosphoBG_5Trans0.75LocFilter.tsv', sep = '\t', header = TRUE)

Lib_combined_peptides <- as.list(unique(Lib_combined$StrippedPeptide))
Lib_12fxn_peptides <- as.list(unique(Lib_12fxn$StrippedPeptide))
Lib_3SS_peptides <- as.list(unique(Lib_3SS$StrippedPeptide))


DIAvLib <- list(
  DIA = threeSS_peptides,
  Lib = Lib_3SS_peptides
  
)


upset(fromList(DIAvLib), order.by = "freq")
ggvenn(DIAvLib) +
  ggtitle('3 DDA Single-Shots')



