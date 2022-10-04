import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

###USER INPUT###

#Uncomment which instrument you're working on data from
platform = 'Pro'
# platform = 'Pro_SmallLibSearch_LocFilter'

report_directory_path = 'Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/directDIA/' + platform + "/"     #Where is your Spectronaut output report?

if platform == 'Pro':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220714_100348_PTMDIAProject_Pro_PhosphoBG_Report.tsv', delimiter = '\t')

if platform == 'Pro_SmallLibSearch_LocFilter':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220902_105134_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_SmallLib0.75Loc_Report.tsv', delimiter = '\t')

#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Lights.tsv', delimiter= '\t')        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
lights['CAPS'] = lights['Peptides'].str.upper()

pattern = r'\[.*?\]'
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Heavies.tsv', delimiter= '\t')['Modified_HeaviesAnnotated'][0:234]
heavies_stripped = [re.sub(pattern, '', s) for s in heavies]


#Finds all pSTY sites in sequence
def get_pos(mod_seq):
    mod_sites = []
    seq = mod_seq.split('_')[1]
    seq = seq.replace('[+57]','').replace('[+16]','')

    pos = seq.find('[+80]')

    while pos != -1:
        mod_sites.append(pos-1)
        seq = seq[:pos] + seq[pos+5:]
        # print(seq)
        pos = seq.find('[+80]')

    return (mod_sites)      #Returns positions of all phosphorylated STY's as index, NOT absolute position in sequence

# print(get_pos('_C[+57]LS[+80]S[+80]IVDSISSEER_'))



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
        # print(len(conc))
        conc = conc.dropna(subset=['EG.PTMAssayProbability'])

        for index, row in conc.iterrows():

            if type(row['EG.PTMProbabilities [Phospho (STY)]']) != float:
                prob = row['EG.PTMProbabilities [Phospho (STY)]'].split(';')
            seq = row['PEP.StrippedSequence']
            mod_seq = row['EG.IntPIMID']

            info = tracker_dict[seq]  # STY site info


            sites = get_pos(mod_seq) # indexed char position

            for s in sites:

                sty_num = info[s]               #Which number sty is it in the sequence

                loc_prob = prob[sty_num-1]
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




# report = localize(heavies_stripped, 'Heavy')
report = localize(lights['CAPS'], 'Lights')

































