import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

wd = 'Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/'

lights = pd.read_csv(wd + 'Modified_Lights.tsv', delimiter= '\t')
heavies = pd.read_csv(wd +  'Modified_Heavies.tsv', delimiter= '\t')


# diann_lights = []
# for index, row in lights.iterrows():
#
#     spectronaut = row['Modified']
#     split = spectronaut.split('_')[1]
#     diann = split.replace('[+80]', '(UniMod:21)')
#     diann_lights.append(diann)
#
#
# lights['DIANN_Mods'] = diann_lights
# lights.to_csv(wd + 'Modified_Lights.tsv', sep= '\t', index = False )

#
# diann_heavies = []
# for index, row in heavies.iterrows():
#
#     spectronaut = str(row['Modified_HeaviesAnnotated'])
#     diann = spectronaut.replace('[+57]','').replace('[+80]','(UniMod:21)').replace('[+8]','(UniMod:259)').replace('[+10]','(UniMod:267)').replace('[+6]','')
#     diann_heavies.append(diann)
#
# heavies['DIANN_Mods'] = diann_heavies
# heavies.to_csv(wd + 'Modified_Heavies.tsv', sep = '\t', index = False)

library_heavies = []
for index, row in heavies.iterrows():
    mod = row['Modified_HeaviesAnnotated']
    if type(mod) == str:
        new = '_'+mod+'_'
    else:
        new = mod
    library_heavies.append(new)

heavies['LibraryHeavies'] = library_heavies
heavies.to_csv(wd+'Modified_Heavies.tsv', sep = '\t', index = False)