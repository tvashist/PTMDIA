require(dplyr)
require(tidyr)
require(stringr)
require(data.table)
library(ggplot2)
library("RColorBrewer")

library <- read.table(file = 'Z:/Helium_Tan/PTMDIAProject_SpectralLibraries/Exports/Pro_3DDASpikedOnly/PTMDIAProject_3DIASpikedPhosphoBG_5Trans0.75LocFilter.tsv', sep = '\t', header = TRUE )
heavies <- read.table(file = 'Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Heavies.tsv', sep = '\t', header = TRUE)



View(library)

annotated <- c()
View(annotated)

for (h in heavies$Modified){
  len = nchar(h)
  new = substring(h, 2, len-1)
  annotated <- append(annotated,new)

  
}


