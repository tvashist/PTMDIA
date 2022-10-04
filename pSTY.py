import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


###USER INPUT###

#Uncomment which instrument you're working on data from
# platform = 'Exploris_FAIMS'
# platform = 'Exploris_NoFAIMS'
platform = 'Pro_SmallLibSearch_LocFilter'
# platform = 'SCP'
# platform = 'Pro_12fxnOnlySearch_LocFilter'

report_directory_path = 'Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/' + platform + "/"     #Where is your Spectronaut output report?

if platform == 'Pro_SmallLibSearch_LocFilter':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220902_105134_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_SmallLib0.75Loc_Report.tsv', on_bad_lines= 'skip', delimiter= '\t',low_memory = False)
    print("READ")


counts = []
for x in spectronaut['PEP.StrippedSequence']:
    s = x.count('S')
    t = x. count('T')
    y = x.count('Y')
    counts.append(s+t+y)

print(mean(counts))