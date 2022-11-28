import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

###USER INPUT###

#Uncomment which instrument you're working on data from
# platform = 'Exploris'
# platform = 'Exploris_FAIMS'
# platform = 'TimsTOF_Pro'
platform = "TimsTOF_SCP"

#Uncommment which library you used
library = 'Library_3SS_Spiked'
# library = 'Library_12fxn_NotSpiked'
# library = 'Library_Combined'
# library = 'directDIA'


report_directory_path = 'Z:/Helium_Tan/FINAL_PTMDIA/' + platform + '/DIANN/' + library + '/SearchOutputs_LocalizationInfo/dia_nn/out/'
diann = pd.read_csv(report_directory_path + 'report.tsv', delimiter='\t', low_memory=False)
print('read')

if platform == 'TimsTOF_SCP':
    conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0]               #Spiked amounts of peptide used for this instrument platform
    norm_spike = 0.2

else:
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    norm_spike = 1.0


run = []
spike = []

for r in diann['Run']:
    parts = r.split('_')

    fmol = [i for i in parts if 'fmol' in i]

    if len(fmol) == 0:      #NoSpikeRuns
        amt = 0

    else:
        fetch = fmol[0].split('fmol')[0]
        amt = fetch.replace('p','.')

    spike.append(float(amt))


    if 'REDO' in r:
        rep = int(parts[-2])

    else:
        rep = int(parts[-1].split('.')[0])

    run.append(rep)
diann['Rep'] = run
diann['Spike'] = spike


#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['DIANN_Mods']


diann_lights = diann.loc[diann['Modified.Sequence'].isin(lights)]         #Only looks at data for the singly-charged non-human light spiked peptides

def localize(spiked_list, type_spike):
    dist = {}
    for x in conditions:

        current = diann_lights[diann_lights['Spike'] == x]
        if x == 10.0:
            current.to_csv(report_directory_path + 'ten.tsv', sep ='\t')
        if x not in dist:
            dist[x] = []
            dist[x] = [n for n in current['PTM.Site.Confidence']]


    df = pd.DataFrame.from_dict(dist, orient='index')
    df = df.transpose()
    df = df.reset_index(drop=True)

    path = report_directory_path + '/' + platform + '_UnfilteredSiteLocalization/'

    # Create the directory if it doesn't already exist
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)

    df.to_csv(path + type_spike + '_UnfilteredSiteLocalization.tsv', sep='\t')
    # print(ten[0])


localize(lights, "Lights")





























