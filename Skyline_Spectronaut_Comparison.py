import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import openpyxl

# pro = pd.read_csv('20220405_152239_PTMDIAProject_Phospho_TimsTOF_Pro_Report.tsv', delimiter= '\t')
#
#
# sub = pro[['R.FileName', 'EG.IntPIMID','EG.ApexRT']]
# sub.to_csv("Pro Filtered List of RTs.tsv", sep = '\t')


skyline = pd.read_excel('PTMDIAProject_Skyline_Heavies_ManualAnnotation.xlsx')
spec_exploris = pd.read_csv('Exploris Filtered List of RTs.tsv', delimiter='\t')
spec_pro = pd.read_csv('Pro Filtered List of RTs.tsv', delimiter= '\t')

print(spec_exploris.columns)
print(spec_pro.columns)
sequences = skyline['Peptide Modified Sequence']
sequences = ['_'+s+'_' for s in sequences]


found_exploris = []
exploris_rt = []
found_pro = []
pro_rt = []

for seq in sequences:
    # print(seq)
    if seq in spec_exploris['EG.IntPIMID'].values:              #If the skyline peptide was found in Exploris
        found_exploris.append('Y')
        select = spec_exploris.loc[spec_exploris['EG.IntPIMID'] == seq]
        rt = select['EG.ApexRT']
        exploris_rt.append(rt)
        # print('found in exploris!')

    if seq not in spec_exploris['EG.IntPIMID'].values:
        found_exploris.append('N')
        exploris_rt.append('N')
        # print('not found in exploris!')

    if seq in spec_pro['EG.IntPIMID'].values:
        found_pro.append('Y')
        select = spec_pro.loc[spec_pro['EG.IntPIMID'] == seq]
        rt = select['EG.ApexRT']
        pro_rt.append(rt)
        # print('found in pro!')

    if seq not in spec_pro['EG.IntPIMID'].values:
        found_pro.append('N')
        pro_rt.append('N')
        # print('not found in pro!')




skyline['Spectronaut_Exploris'] = found_exploris
skyline['Spectronaut_Pro'] = found_pro
# skyline['Spectronaut_Exploris_RT'] = exploris_rt
# skyline['Spectronaut_Pro_RT'] = pro_rt

skyline.to_csv('Comparing Spectronaut and Skyline.tsv', sep = '\t')


