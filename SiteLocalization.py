import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

###USER INPUT###

#Uncomment which instrument you're working on data from
platform = 'Pro'

report_directory_path = 'Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/directDIA/' + platform + "/"     #Where is your Spectronaut output report?

if platform == 'Pro':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220714_100348_PTMDIAProject_Pro_PhosphoBG_Report.tsv', delimiter = '\t')
    print('READ')

#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Lights.tsv', delimiter= '\t')        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
lights['CAPS'] = lights['Peptides'].str.upper()

pattern = r'\[.*?\]'
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Heavies.tsv', delimiter= '\t')['Modified_HeaviesAnnotated'][0:234]
heavies_stripped = [re.sub(pattern, '', s) for s in heavies]


def get_pos(mod_seq):
    seq = mod_seq.split('_')[1]
    mod_site = seq.find('[')
    return int(mod_site)

def localize(spiked_list, type_spike):

    tracker_dict = {}

    for seq in spiked_list:

        tracker = [None for i in range(0, len(seq))]
        num = 0

        for i in range(0, len(seq)):
            if seq[i] == "S" or seq[i] == "T" or seq[i] == "Y":
                num += 1
                tracker[i] = num

            else:
                tracker[i] = 0

        if seq not in tracker_dict:
            tracker_dict[seq] = None
            tracker_dict[seq] = tracker

        else:
            tracker_dict[seq] = tracker



    report = {}
    for c in conditions:

        probabilities = []

        spec_conc = spectronaut.loc[spectronaut['R.Condition'].str.contains(str(c))]
        if type_spike == 'Heavy':
            spec_conc = spec_conc.loc[spec_conc['EG.PTMLocalizationProbabilities'].str.contains(str('Label'))]      #For heavy spikes, look only at heavy versions of the peptides
        conc = spec_conc.loc[spec_conc['PEP.StrippedSequence'].isin(spiked_list)]
        conc = conc.dropna(subset=['EG.PTMAssayProbability'])

        for index, row in conc.iterrows():

            if type(row['EG.PTMProbabilities [Phospho (STY)]']) != float:
                prob = row['EG.PTMProbabilities [Phospho (STY)]'].split(';')
            seq = row['PEP.StrippedSequence']
            mod_seq = row['EG.IntPIMID']

            info = tracker_dict[seq]  # STY site info

            site = get_pos(mod_seq) - 1  # indexed char position
            # print(site)

            sty_num = info[site] - 1

            loc_prob = prob[sty_num]
            probabilities.append(loc_prob)

        k = str(c) + ' fmol'
        if k not in report:
            report[k] = None
            report[k] = probabilities
        else:
            report[k] = probabilities

    # print(dict)
    df = pd.DataFrame.from_dict(report, orient='index')
    df = df.transpose()
    df = df.reset_index(drop=True)

    df.to_csv(report_directory_path + type_spike + '_SiteLocalization.tsv', sep='\t')




report = localize(heavies_stripped, 'Heavy')


































