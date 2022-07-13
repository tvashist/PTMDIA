import pandas
import pandas as pd
import re
import math
from statistics import mean,stdev,mode

#Select which instrument you're working on data from
# platform = 'Exploris'
# platform = 'Pro'
platform = 'SCP'

if platform == 'SCP':
    conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0]
    spectronaut = pd.read_csv('20220618_114730_PTMDIAProject_SCP_PhosphoBGCurve_Report.tsv', delimiter='\t', low_memory= False)


if platform == 'Pro':
    conditions = [0.04, 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv('20220405_152239_PTMDIAProject_Phospho_TimsTOF_Pro_Report.tsv', delimiter= '\t', low_memory= False)


if platform == 'Exploris':
    conditions = [0.04, 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv('20220411_084730_PTMDIAProject_PhoshpoTestCurve_DIA_Exploris_Report.tsv', delimiter = '\t', low_memory= False)

#
# ###Lights###
# summary_lights = {}
# lights = pd.read_csv('Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
#
#
# for sequence in lights:
#
#
#     find = sequence     #Spiked in peptide sequence
#
#
#     single = spectronaut.loc[spectronaut['EG.IntPIMID'] == find]        #Only look at the entries where that specific peptide was found
#     single = single[['R.Condition','R.Replicate','EG.IntPIMID','FG.Quantity']]      #Keep only the column indicating condition/spike, peptide sequence, and quantity
#     single = single[single['FG.Quantity'].notnull()]                    #Remove columns where quantity is null
#
#     #Add new column with spiked amount
#     spiked = []
#     for x in single['R.Condition']:
#         amount = float(x.replace('fmol',''))
#         spiked.append(amount)
#
#     single['Spiked'] = spiked
#
#     # print(single)
#
#
#
#     df_list = []
#     for x in conditions:                                    #For every point on the curve, create a new data frame
#
#         spike_label = str(x) + 'fmol'
#         if spike_label not in summary_lights:
#             summary_lights[spike_label] = {}
#             summary_lights[spike_label]['Reps'] = []
#             summary_lights[spike_label]['CVs'] = []
#
#         data = {}
#         current = single[single['Spiked'] == x]             #New df where you are only looking at entries associated with that spiked amount
#
#         data['Type Phosphopeptide'] = 'Light Non-Human'     #Annotate whether we spiked in non-human light or human heavy
#         data['Peptide'] = find
#         data['Spike'] = x
#         data['Log Spike'] = math.log10(x)                   #Plot log spiked amount vs log quantity
#
#         #If spiked peptide is not found at this point on the curve
#         if len(current) == 0:
#             data['Number of Reps'] = 0
#             data['Mean Quant'] = None
#             data['Log Mean Quant'] = None
#             data['Percent CV'] = None
#             data['Expected Ratio'] = conditions[-4] / x
#
#
#         #If spiked peptide is found at this point on the curve
#         if len(current) > 0:
#             grouped = current.groupby('R.Replicate', as_index= False).agg({'R.Condition':'first','Spiked':'first','EG.IntPIMID':'first','FG.Quantity':'sum' })      #Condense to one line per replicate, sum quantities
#             num_reps = len(grouped['R.Replicate'])              #Number of replicates it's found in
#             quant = math.log10(mean(grouped['FG.Quantity']))    #Take the log of the mean quantity
#             cv = (grouped['FG.Quantity'].std() / mean(grouped['FG.Quantity'])) * 100        #Calculate percent CV between quantities found across replicates
#
#             #Add relevant data to dictionary
#             data['Quant_EachRep'] = [x for x in grouped['FG.Quantity']]
#             data['Number of Reps'] = num_reps
#             data['Mean Quant'] = mean(grouped['FG.Quantity'])
#             data['Log Mean Quant'] = quant
#             data['Percent CV'] = cv
#             data['Expected Ratio'] = conditions[-4] / x
#
#             #Add reps and CVs information to summarizing table
#             summary_lights[spike_label]['Reps'].append(num_reps)
#             summary_lights[spike_label]['CVs'].append(cv)
#
#
#             df_list.append(data)  # Collect all dictionaries into a list, so by the end it will be a list of dicts for each point on the curve
#
#         df = pd.DataFrame(df_list)  # Combine all dicts into one data frame and export separately for each spiked peptide
#
#
# ###RATIOS### This is done where each injection a peptide is found in, the quant is ratio'd to the mean quant at the 'set' point
#     if len(df) != 0:
#         ratios = []
#         set = df[df['Spike'] == 0.2]
#         print(set)
#         set_mean = mean(set['Quant_EachRep'].tolist()[0])
#
#         for index, row in df.iterrows():
#             current = row['Quant_EachRep']
#             if current != None:
#                 rat = [set_mean/x for x in current]
#                 ratios.append(rat)
#             else:
#                 ratios.append(None)
#
#         df['Actual Ratios'] = ratios
#
#
#         new_df = pd.DataFrame(columns = ['Peptide', 'Spike', 'Expected Ratio', 'Actual Ratio'])
#         for index, row in df.iterrows():
#             peptide = row['Peptide']
#             spike = row['Spike']
#             expected = row['Expected Ratio']
#
#             for x in row['Actual Ratios']:
#                 sub = {'Peptide': peptide, 'Spike': spike, 'Expected Ratio': expected, 'Actual Ratio': x}
#                 sub = pd.DataFrame.from_dict([sub])
#                 new_df = pandas.concat([new_df,sub])
#
#
#         new_df.to_csv('PhosphoBG_Curve/SCP_ManiPlotting/' + platform + '_Lights_outputs_Found/' + find+'_replicates_output.tsv', sep = '\t')

###Heavies###
#All the same comments applied to the lights code above apply to heavies below
summary_heavies = {}

heavies = pd.read_csv('Modified_Heavies.tsv', delimiter= '\t')['Modified']
found_heavies = spectronaut.loc[spectronaut['EG.PTMLocalizationProbabilities'].str.contains('Label')]        #Only look at the entries that contain heavy modification


for sequence in heavies:
    found_at_spike = 0
    PeptideFound = None        #boolean tracking if peptide was found or not
    find = sequence     #Spiked in peptide sequence
    # print(find)


    single = found_heavies.loc[found_heavies['EG.IntPIMID'] == find]        #Only look at the entries where that specific peptide was found
    single = single[['R.Condition','R.Replicate','EG.IntPIMID','FG.Quantity']]
    single = single[single['FG.Quantity'].notnull()]



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
        set = df[df['Spike'] == 0.2]
        if int(set['Number of Reps']) > 0:               #If there is quant at this spike level
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


            new_df.to_csv('PhosphoBG_Curve/SCP_ManiPlotting/' + platform + '_Heavies_outputs_Found/' + find+'_replicates_output.tsv', sep = '\t')



















