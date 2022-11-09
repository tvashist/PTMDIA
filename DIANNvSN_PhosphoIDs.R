require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library(readr)
library(UpSetR)

#Load Pro 12fxn search from DIANN and Spectronaut

diann_12fxn <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/DIANN/Library_12fxn_NotSpiked/SearchOutputs/dia_nn/out/report.tsv")
spectronaut_12fxn <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_12fxn_NotSpiked/SearchOutputs/20220906_093527_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_12fxnLib0.75Loc_Report.tsv")

diann_combined <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/DIANN/Library_Combined/SearchOutputs/dia_nn/out/report.tsv")
spectronaut_combined <- read_tsv("Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_Combined/SearchOutputs/20220803_134456_PTMDIAProject_DIACurveAnalysis_WithSpecLib_Report_addedFGLabel.tsv")

#Select only "no spike" runs
diann_nospike <- diann_12fxn %>% filter(grepl("NoSpike",Run))
spec_nospike <- spectronaut_12fxn %>% filter(grepl("NoSpike",R.FileName))


diann_modified <- diann_nospike$Modified.Sequence %>% gsub("UniMod:1","", .) %>% gsub("UniMod:4", "", .) %>% gsub("UniMod:35","", .) %>% gsub("UniMod:21","[+80]", .) %>% gsub("[()]", "", .)
spec_modified <- spec_nospike$EG.IntPIMID %>% gsub("[+57]", "", ., fixed=TRUE) %>% gsub("[+42]","",.,fixed = TRUE) %>% gsub("[+16]", "", ., fixed = TRUE) %>% gsub("_","",.) 



make_upset <- function(diann, spectronaut) {
  
  diann_sequences <- unique(diann)
  spectronaut_sequences <- unique(spectronaut)
  
  
  overlap <- list(DIANN = diann_sequences, SN = spectronaut_sequences)
  
  u <- upset(fromList(overlap),sets.bar.color=c("maroon","blue"))

  
}


p <- make_upset(diann_modified, spec_modified)
p



                                                                                                  



