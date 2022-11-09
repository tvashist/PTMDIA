import pandas
import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

def makePath(filepath):
    isExist = os.path.exists(filepath)
    if not isExist:
        os.makedirs(filepath)


#Uncomment which instrument you're working on data from#
# platform = 'Exploris'
# platform = 'Exploris_FAIMS'
platform = 'TimsTOF_Pro'
# platform = 'TimsTOF_SCP'

#Library#
# library = 'Library_3SS_Spiked'
# library = 'Library_12fxn_NotSpiked'
# library = 'Library_Combined'
library = 'directDIA'

if platform == 'TimsTOF_SCP':
    conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0]               #Spiked amounts of peptide used for this instrument platform
    norm_spike = 0.2

else:
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    norm_spike = 1.0


report_directory_path = 'Z:/Helium_Tan/FINAL_PTMDIA/' + platform + '/Spectronaut/' + library + '/SearchOutputs/'
report = '20221027_100139_PTMDIAProject_TimsTOFPro_directDIA_DIACurveAnalysis_Report.tsv'
spectronaut = pd.read_csv(report_directory_path + report, delimiter='\t', low_memory=False)
print('read')



# ###Lights###
summary_lights = {}
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut


for sequence in lights:


    find = sequence     #Spiked in peptide sequence


    single = spectronaut.loc[spectronaut['EG.IntPIMID'] == find]        #Only look at the entries where that specific peptide was found
    single = single[['R.Condition','R.Replicate','EG.IntPIMID','FG.Quantity']]      #Keep only the column indicating condition/spike, peptide sequence, and quantity
    single = single[single['FG.Quantity'].notnull()]                    #Remove columns where quantity is null

    #Add new column with spiked amount
    spiked = []
    for x in single['R.Condition']:
        amount = float(x.replace('fmol',''))
        spiked.append(amount)

    single['Spiked'] = spiked

    # print(single)



    df_list = []
    for x in conditions:                                    #For every point on the curve, create a new data frame

        spike_label = str(x) + 'fmol'
        if spike_label not in summary_lights:
            summary_lights[spike_label] = {}
            summary_lights[spike_label]['Reps'] = []
            summary_lights[spike_label]['CVs'] = []

        data = {}
        current = single[single['Spiked'] == x]             #New df where you are only looking at entries associated with that spiked amount

        data['Type Phosphopeptide'] = 'Light Non-Human'     #Annotate whether we spiked in non-human light or human heavy
        data['Peptide'] = find
        data['Spike'] = x
        data['Log Spike'] = math.log10(x)                   #Plot log spiked amount vs log quantity

        #If spiked peptide is not found at this point on the curve
        if len(current) == 0:
            data['Quant_EachRep'] = []
            data['Number of Reps'] = 0
            data['Mean Quant'] = None
            data['Log Mean Quant'] = None
            data['Percent CV'] = None
            data['Expected Ratio'] = conditions[-4] / x


        #If spiked peptide is found at this point on the curve
        if len(current) > 0:
            grouped = current.groupby('R.Replicate', as_index= False).agg({'R.Condition':'first','Spiked':'first','EG.IntPIMID':'first','FG.Quantity':'sum' })      #Condense to one line per replicate, sum quantities
            num_reps = len(grouped['R.Replicate'])              #Number of replicates it's found in
            quant = math.log10(mean(grouped['FG.Quantity']))    #Take the log of the mean quantity
            cv = (grouped['FG.Quantity'].std() / mean(grouped['FG.Quantity'])) * 100        #Calculate percent CV between quantities found across replicates

            #Add relevant data to dictionary
            data['Quant_EachRep'] = [x for x in grouped['FG.Quantity']]
            data['Number of Reps'] = num_reps
            data['Mean Quant'] = mean(grouped['FG.Quantity'])
            data['Log Mean Quant'] = quant
            data['Percent CV'] = cv
            data['Expected Ratio'] = conditions[-4] / x

            #Add reps and CVs information to summarizing table
            summary_lights[spike_label]['Reps'].append(num_reps)
            summary_lights[spike_label]['CVs'].append(cv)


            df_list.append(data)  # Collect all dictionaries into a list, so by the end it will be a list of dicts for each point on the curve

        df = pd.DataFrame(df_list)  # Combine all dicts into one data frame and export separately for each spiked peptide


###RATIOS### This is done where each injection a peptide is found in, the quant is ratio'd to the mean quant at the 'set' point
    if len(df) != 0:
        ratios = []
        set = df[df['Spike'] == norm_spike]
        # print(set)
        if len(set) > 0:               #If there is quant at this spike level
            set_mean = mean(set['Quant_EachRep'].tolist()[0])

            for index, row in df.iterrows():
                current = row['Quant_EachRep']
                if isinstance(current, float) == False:
                    rat = [set_mean/x for x in current]
                    ratios.append(rat)
                else:
                    ratios.append(None)

            df['Actual Ratios'] = ratios


            new_df = pd.DataFrame(columns = ['Peptide', 'Spike', 'Expected Ratio', 'Actual Ratio'])
            for index, row in df.iterrows():
                peptide = row['Peptide']
                spike = row['Spike']
                expected = row['Expected Ratio']

                if row['Actual Ratios'] != None:
                    for x in row['Actual Ratios']:
                        sub = {'Peptide': peptide, 'Spike': spike, 'Expected Ratio': expected, 'Actual Ratio': x}
                        sub = pd.DataFrame.from_dict([sub])
                        new_df = pandas.concat([new_df,sub])



            path = report_directory_path + platform + '_Lights_outputs_Found/PlottingReplicates/'
            all_spikes_path = report_directory_path + platform + '_AllSpikes_Replicates/'

            makePath(path)
            makePath(all_spikes_path)


            new_df.to_csv(path + find + '_replicates_output.tsv', sep='\t')
            new_df.to_csv(all_spikes_path + find + '_replicates_output.tsv', sep='\t')






###Heavies###
#All the same comments applied to the lights code above apply to heavies below
summary_heavies = {}

heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['Modified_HeaviesAnnotated'][0:234]
pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv',delimiter='\t')['Modified_HeaviesAnnotated'][0:234]
found_heavies = spectronaut.loc[spectronaut['FG.LabeledSequence'].str.contains('Label')]  # Only look at the entries that contain heavy modification

new_annotations = []
for seq in found_heavies['FG.LabeledSequence']:  # This column contains sequences with all mod annotations including heavy labeling
    new = seq.replace('Phospho (STY)', '+80').replace('Carbamidomethyl (C)', '+57').replace('Label:13C(6)15N(4)','+10').replace(
        'Label:13C(6)15N(2)', '+8').replace('[Oxidation (M)]', '').replace('[Label:13C(5)15N(1)]', '').replace(
        '[Label:13C(6)15N(2)]', '').replace('[Label:T-5]', '').replace('[Label:13C(6)15N(1)]', '').replace(
        '[Acetyl (Protein N-term)]', '')
    new = new[1:-1]
    new_annotations.append(new)
found_heavies['AnnotatedSequence'] = new_annotations

for sequence in heavies:
    find = sequence  # Spiked in peptide sequence
    found_at_spike = 0
    PeptideFound = None

    single = found_heavies.loc[found_heavies['AnnotatedSequence'] == find]  # Only look at the entries where that specific peptide was found
    single = single[['R.Condition', 'R.Replicate', 'EG.IntPIMID', 'FG.Quantity', 'AnnotatedSequence']]
    single = single[single['FG.Quantity'].notnull()]
# summary_heavies = {}
#
# heavies = pd.read_csv('Z:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Heavies.tsv', delimiter= '\t')['Modified']
# found_heavies = spectronaut.loc[spectronaut['EG.PTMLocalizationProbabilities'].str.contains('Label')]        #Only look at the entries that contain heavy modification
#
#
# for sequence in heavies:
#     found_at_spike = 0
#     PeptideFound = None        #boolean tracking if peptide was found or not
#     find = sequence     #Spiked in peptide sequence
#
#
#     single = found_heavies.loc[found_heavies['EG.IntPIMID'] == find]        #Only look at the entries where that specific peptide was found
#     single = single[['R.Condition','R.Replicate','EG.IntPIMID','FG.Quantity']]
#     single = single[single['FG.Quantity'].notnull()]



    spiked = []
    for x in single['R.Condition']:
        amount = float(x.replace('fmol',''))
        spiked.append(amount)

    single['Spiked'] = spiked

    df_list = []
    for x in conditions:

        spike_label = str(x) + 'fmol'
        if spike_label not in summary_heavies:
            summary_heavies[spike_label] = {}
            summary_heavies[spike_label]['Reps'] = []
            summary_heavies[spike_label]['CVs'] = []



        data = {}
        current = single[single['Spiked'] == x]           # new df where you are only looking at entries associated with that spiked amount

        data['Type Phosphopeptide'] = 'Heavy Human'
        data['Peptide'] = find
        data['Spike'] = x
        data['Log Spike'] = math.log10(x)

        if len(current) == 0:
            data['Quant_EachRep'] = []
            data['Number of Reps'] = 0
            data['Mean Quant'] = None
            data['Log Mean Quant'] = None
            data['Percent CV'] = None
            data['Expected Ratio'] = conditions[-4] / x



        if len(current) > 0:
            found_at_spike += 1
            grouped = current.groupby('R.Replicate', as_index= False).agg({'R.Condition':'first','Spiked':'first','EG.IntPIMID':'first','FG.Quantity':'sum' })
            num_reps = len(grouped['R.Replicate'])
            quant = math.log10(mean(grouped['FG.Quantity']))
            cv = (grouped['FG.Quantity'].std() / mean(grouped['FG.Quantity'])) * 100

            data['Quant_EachRep'] = [x for x in grouped['FG.Quantity']]
            data['Number of Reps'] = num_reps
            data['Mean Quant'] = mean(grouped['FG.Quantity'])
            data['Log Mean Quant'] = quant
            data['Percent CV'] = cv
            data['Expected Ratio'] = conditions[-4] / x
            data['Type Phosphopeptide'] = 'Heavy Human'


            #Add reps and CVs information to summarizing table
            summary_heavies[spike_label]['Reps'].append(num_reps)
            summary_heavies[spike_label]['CVs'].append(cv)

        df_list.append(data)

    df = pd.DataFrame(df_list)
    if found_at_spike > 0:      #If the peptide is detected at any spike level, this counter should be > 0. Set boolean to true.
        Found = True
    else:
        Found = False


###HEAVY RATIOS### This is done where each injection a peptide is found in, the quant is ratio'd to the mean quant at the 'set' point
    if Found == True:
        ratios = []
        set = df[df['Spike'] == norm_spike]

        quants = set['Quant_EachRep'].tolist()[0]
        if len(quants) > 0:               #If there is quant at this spike level
            set_mean = mean(set['Quant_EachRep'].tolist()[0])

            for index, row in df.iterrows():
                current = row['Quant_EachRep']
                if isinstance(current, float) == False:
                    rat = [set_mean/x for x in current]
                    ratios.append(rat)
                else:
                    ratios.append(None)

            df['Actual Ratios'] = ratios


            new_df = pd.DataFrame(columns = ['Peptide', 'Spike', 'Expected Ratio', 'Actual Ratio'])
            for index, row in df.iterrows():
                peptide = row['Peptide']
                spike = row['Spike']
                expected = row['Expected Ratio']

                if row['Actual Ratios'] != None:
                    for x in row['Actual Ratios']:
                        sub = {'Peptide': peptide, 'Spike': spike, 'Expected Ratio': expected, 'Actual Ratio': x}
                        sub = pd.DataFrame.from_dict([sub])
                        new_df = pandas.concat([new_df,sub])

            path = report_directory_path + platform + '_Heavies_outputs_Found/PlottingReplicates/'
            all_spikes_path = report_directory_path + platform + '_AllSpikes_Replicates/'

            makePath(path)
            makePath(all_spikes_path)

            new_df.to_csv(path + find + '_replicates_output.tsv', sep='\t')
            new_df.to_csv(all_spikes_path + find + '_replicates_output.tsv', sep='\t')






































