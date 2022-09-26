import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

###USER INPUT###

#Uncomment which instrument you're working on data from
# platform = 'Exploris_FAIMS'
platform = 'Exploris_NoFAIMS'
# platform = 'Pro_SmallLibSearch_LocFilter'
# platform = 'SCP'
# platform = 'Pro_12fxnOnlySearch_LocFilter'

report_directory_path = 'Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/directDIA/' + platform + "/"     #Where is your Spectronaut output report?

if platform == 'Exploris_NoFAIMS':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220921_094453_PTMDIAProject_Exploris_DIACurveAnalysis_directDIA_Report.tsv', delimiter = '\t')

#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut

for sequence in lights[0:1]:


    find = sequence     #Spiked in peptide sequence


    single = spectronaut.loc[spectronaut['EG.IntPIMID'] == find]            #Only look at the entries where that specific peptide was found
    single.to_csv(report_directory_path + 'single.tsv', sep = '\t')