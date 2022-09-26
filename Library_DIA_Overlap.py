import math
from statistics import mean,stdev,mode
import os
import pandas as pd

library = pd.read_csv('Z:/Helium_Tan/PTMDIAProject_SpectralLibraries/Pro_12fxnOnly/PTMDIAProject_TimsTOFPro_12fxnOnly.tsv', delimiter= '\t',low_memory = False)
 # print(len(library))
phospho_library = library[library['IntLabeledPeptide'].str.contains('80')]
# print(len(phospho_library))

DIACurve = pd.read_csv('Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro_12fxnOnlySearch_LocFilter/20220906_093527_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_12fxnLib0.75Loc_Report.tsv', delimiter= '\t',low_memory = False)
phospho_DIACurve = DIACurve.loc[DIACurve['EG.IntPIMID'].str.contains('80')]

library_peps = set(list(phospho_library['StrippedPeptide']))
DIA_peps = set(list(phospho_DIACurve['PEP.StrippedSequence']))
overlap = library_peps.intersection(DIA_peps)

print('Library contains: ' + str(len(library_peps)) + " peptides")
print('DIA curve contains: ' + str(len(DIA_peps)) + " peptides")
print("There are " + str(len(overlap)) + " peptides in both sets")

DIA_only = DIA_peps - library_peps

scores = {}
for index, row in phospho_DIACurve.iterrows():
    seq = row['PEP.StrippedSequence']
    cscore = row['EG.Cscore']

    if seq in DIA_only:
        if seq not in scores:
            scores[seq] = []
            scores[seq].append(cscore)
        else:
            scores[seq].append(cscore)

print(scores)




