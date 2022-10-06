import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

wd = 'Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/'

lights = pd.read_csv(wd + 'Modified_Lights.tsv', delimiter= '\t')


# site = []
# for index, row in lights.iterrows():
#     mod = row['Peptides']
#     for i in range(0,len(mod)):
#         print(mod[i])
#         if mod[i].islower():
#             site.append(i+1)
#
# lights['Mod Site'] = site
# lights.to_csv(wd +"Modified_Lights_SitePositions.tsv", sep = '\t')

diann_mods = []
for index, row in lights.iterrows():

    spectronaut = row['Modified']
    split = spectronaut.split('_')[1]
    diann = split.replace('[+80]', '(UniMod:21)')
    diann_mods.append(diann)


lights['DIANN_Mods'] = diann_mods
lights.to_csv(wd + 'Modified_Lights.tsv', sep= '\t', index = False )

