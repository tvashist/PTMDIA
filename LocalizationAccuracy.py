import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

###USER INPUT###

#Uncomment which instrument you're working on data from
# platform = 'Exploris'
platform = 'Exploris_FAIMS'
# platform = 'TimsTOF_Pro'
# platform = "TimsTOF_SCP"

#Uncommment which library you used
library = 'Library_3SS_Spiked'
# library = 'Library_12fxn_NotSpiked'
# library = 'Library_Combined'
# library = 'directDIA'


report_directory_path = "Z:/Helium_Tan/FINAL_PTMDIA/" + platform + "/Spectronaut/" + library + "/SearchOutputs/"     #Where is your Spectronaut output report?
report = '20221101_145814_PTMDIAProject_ExplorisFAIMS_3SSLib_DIACurveAnalysis_Report.tsv'
spectronaut = pd.read_csv(report_directory_path + report, sep=  '\t', low_memory= False )
print('read')


if platform == 'TimsTOF_SCP':
    conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0]

else:
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]


#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['Modified']



def localize(spiked_list, type_spike):

    dist = {}
    localized = {}


    for c in conditions:

        probabilities = []

        spec_conc = spectronaut.loc[spectronaut['R.Condition'].str.contains(str(c))]
        if type_spike == 'Heavy':
            spec_conc = spec_conc.loc[spec_conc['EG.PTMLocalizationProbabilities'].str.contains(str('Label'))]      #For heavy spikes, look only at heavy versions of the peptides
        conc = spec_conc.loc[spec_conc['EG.IntPIMID'].isin(spiked_list)]
        # print(len(conc))
        conc = conc.dropna(subset=['EG.PTMAssayProbability'])

        for index, row in conc.iterrows():

            if type(row['EG.PTMProbabilities [Phospho (STY)]']) != float:
                prob = row['EG.PTMProbabilities [Phospho (STY)]'].split(';')

            seq = row['PEP.StrippedSequence']
            mod_seq = row['EG.IntPIMID']

            #Finding probabilitiy info

            all_sty = []
            find_sty = re.finditer('[STY](\[\+80\])?', mod_seq)
            for match in find_sty:
                all_sty.append(match.start())

            all_phos = []
            find_phos = re.finditer('[STY](\[\+80\])', mod_seq)
            for match in find_phos:
                all_phos.append(match.start())

            if len(all_sty) == len(prob):
                for x in range(0, len(all_sty)):
                    if all_sty[x] in all_phos:
                        probabilities.append(prob[x])

        loc = 0
        total = len(probabilities)
        for x in probabilities:
            if float(x) >= 0.75:
                loc += 1

        perc_loc = (loc/total) * 100
        print(perc_loc)




        k = str(c) + ' fmol'
        if k not in report:
            dist[k] = None
            dist[k] = probabilities

        else:
            dist[k] = probabilities


        if k not in localized:
            localized[k] = None
            localized[k] = perc_loc

        else:
            localized[k] = perc_loc


    print(localized)
    # print(dict)
    df = pd.DataFrame.from_dict(dist, orient='index')
    df = df.transpose()
    df = df.reset_index(drop=True)

    df_loc = pd.DataFrame.from_dict(localized, orient = 'index')
    df_loc = df_loc.transpose()
    df_loc = df_loc.reset_index(drop = True)



    path = report_directory_path + '/' + platform + '_UnfilteredSiteLocalization/'

    # Create the directory if it doesn't already exist
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)

    df.to_csv(path + type_spike + '_UnfilteredSiteLocalization.tsv', sep='\t')
    df_loc.to_csv(path + type_spike + '_PercLocalized.tsv', sep = '\t')






heavies_report = localize(heavies, 'Heavy')
lights_report = localize(lights, 'Lights')































