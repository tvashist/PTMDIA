import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

wd = "Z:/Helium_Tan/FINAL_PTMDIA/TimsTOF_Pro/Spectronaut/Library_Combined/LibraryExport/"
combined_library = pd.read_csv(wd + "PTMDIAProject_TimsTOFPro_Combined_Library.tsv", delimiter='\t')

reference = []
for index, row in combined_library.iterrows():
    ref = row['ReferenceRun']
    if 'SpectralLib' in ref:
        reference.append('12fxn')
    if '50ugPhosphoBG' in ref:
        reference.append('3SS')

combined_library['Source'] = reference

combined_library.to_csv(wd + "Annotated_PTMDIAProject_TimsTOFPro_Combined_Library.tsv", sep= '\t')
