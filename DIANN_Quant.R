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

setwd('Z:/Helium_Tan/DIA_NN/Pro_3DDALib/dia_nn/out/')

lights<- read.table('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Lights.tsv', sep ='\t', header = TRUE, quote = '')
lights <- lights %>% mutate('CAPS' = toupper(Peptides))
View(lights$CAPS)
View(lights$Modified)


report <- read_tsv('report.tsv')
View(head(report))
filtered <- report %>% filter(grepl('UniMod:21', Modified.Sequence)) %>% filter(Stripped.Sequence %in% lights$CAPS) 
filtered$Modified.Sequence <- gsub('UniMod:21', '+80', filtered$Modified.Sequence) %>% gsub("[(UniMod:21)]", "boo", filtered$Modified.Sequence) #Need to replace the parentheses with square brackets
View(filtered)


# 
# lib <- read_tsv('spect_lib.tsv')
# filtered_lib <- lib %>% filter(grepl('UniMod:21', FullUniModPeptideName))
# View(filtered_lib)

# 
# lib <- read_tsv('Z:/Helium_Tan/PTMDIAProject_SpectralLibraries/Exports/Pro_3DDASpikedOnly/PTMDIAProject_3DIASpikedPhosphoBG_5Trans0.75LocFilter.tsv')
# View(head(lib))
# 
# heavy <- lib %>% filter(grepl("(\\[.+8\\])", IntLabeledPeptide))
# View(heavy)

