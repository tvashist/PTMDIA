import pandas as pd

print('hello')
file = pd.read_csv('Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/SCP/20220813_024041_PTMDIAProject_SCP_DIACurveAnalysis_WithSpecLib_Report.tsv', delimiter = '\t', low_memory = False)
print("READ")
short = file.loc[file['R.FileName'].str.contains('2.0fmol')]

short.to_csv('Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/SCP/2fmolshort.tsv', sep = '\t')