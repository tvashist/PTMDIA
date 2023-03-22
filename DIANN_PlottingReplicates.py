#!/usr/bin/env python
# coding: utf-8

# In[2]:


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[3]:


import pandas
import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


# In[4]:


def makePath(filepath):
    isExist = os.path.exists(filepath)
    if not isExist:
        os.makedirs(filepath)


# In[64]:


report_directory_path = 'Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Predicted_Library\dia_nn\out\\'
diann = pd.read_csv(report_directory_path + 'report.tsv', delimiter='\t', low_memory=False)

print('read')


# In[66]:


check = diann.loc[diann['Modified.Sequence'].str.contains('259')]
print(len(check))


# In[68]:


platform = 'Pro'
workflow = 'diann_Predicted'
spikes = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
# spikes = [str(x).replace('.','p') for x in spikes]

if platform == 'SCP':
    conditions = spikes[:9]
else:
    conditions = spikes[4:]
    
norm_spike = 1.0
    
print(conditions)


# In[69]:


lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['DIANN_Mods']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['DIANN_Mods']


# In[71]:


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


# In[72]:


summary_lights = {}


for sequence in lights:


    find = str(sequence)     #Spiked in peptide sequence


    single = diann.loc[diann['Modified.Sequence'] == find]                     #Only look at the entries where that specific peptide was found
    single = single[['Rep','Spike','Modified.Sequence','Precursor.Quantity']]  #Keep only the column indicating condition/spike, peptide sequence, and quantity
    single = single[single['Precursor.Quantity'].notnull()]                    #Remove columns where quantity is null
    single = single[single['Precursor.Quantity'] > 0]                          #Only keep columns where quantity > 0


    df_list = []
    for x in conditions:                                    #For every point on the curve, create a new data frame

        spike_label = str(x) + 'fmol'
        if spike_label not in summary_lights:
            summary_lights[spike_label] = {}
            summary_lights[spike_label]['Reps'] = []
            summary_lights[spike_label]['CVs'] = []

        data = {}
        current = single[single['Spike'] == x]             #New df where you are only looking at entries associated with that spiked amount
         
        data['Workflow'] = workflow
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
            grouped = current.groupby('Rep', as_index= False).agg({'Spike':'first','Modified.Sequence':'first','Precursor.Quantity':'sum' })       #Condense to one line per replicate, sum quantities
            num_reps = len(grouped['Precursor.Quantity'])              #Number of replicates it's found in
            quant = math.log10(mean(grouped['Precursor.Quantity']))    #Take the log of the mean quantity
            cv = (grouped['Precursor.Quantity'].std() / mean(grouped['Precursor.Quantity'])) * 100        #Calculate percent CV between raw quantities found across replicates

            #Add relevant data to dictionary
            data['Quant_EachRep'] = [x for x in grouped['Precursor.Quantity']]
            data['Number of Reps'] = num_reps
            data['Mean Quant'] = mean(grouped['Precursor.Quantity'])
            data['Log Mean Quant'] = quant
            data['Percent CV'] = cv
            data['Expected Ratio'] = conditions[-4] / x

            #Add reps and CVs information to summarizing table
            summary_lights[spike_label]['Reps'].append(num_reps)
            summary_lights[spike_label]['CVs'].append(cv)


            df_list.append(data)  # Collect all dictionaries into a list, so by the end it will be a list of dicts for each point on the curve

    df = pd.DataFrame(df_list)  # Combine all dicts into one data frame and export separately for each spiked peptide
        
    info_path = report_directory_path + platform + '_AllSpikes_Info/'
    makePath(info_path)
    x = find.replace(':', '')
    df.to_csv(info_path + x + '_' + workflow + '_output.tsv', sep='\t')


###RATIOS### This is done where each injection a peptide is found in, the quant is ratio'd to the mean quant at the 'set' point
    if len(df) != 0:
        ratios = []
        set = df[df['Spike'] == norm_spike]

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

            x = find.replace(':', '')
            new_df.to_csv(path + x + '_replicates_output.tsv', sep='\t')
            new_df.to_csv(all_spikes_path + x + '_replicates_output.tsv', sep='\t')


# In[73]:


###Heavies###
#All the same comments applied to the lights code above apply to heavies below
summary_heavies = {}


for sequence in heavies:
    find = sequence  # Spiked in peptide sequence
    found_at_spike = 0
    PeptideFound = None

    single = diann.loc[diann['Modified.Sequence'] == find]            #Only look at the entries where that specific peptide was found
    single = single[['Rep','Spike','Modified.Sequence','Precursor.Quantity']]      #Keep only the column indicating condition/spike, peptide sequence, and quantity
    single = single[single['Precursor.Quantity'].notnull()]                        #Remove columns where quantity is null
    single = single[single['Precursor.Quantity'] > 0]                          #Only keep columns where quantity > 0


    df_list = []
    for x in conditions:

        spike_label = str(x) + 'fmol'
        if spike_label not in summary_heavies:
            summary_heavies[spike_label] = {}
            summary_heavies[spike_label]['Reps'] = []
            summary_heavies[spike_label]['CVs'] = []



        data = {}
        current = single[single['Spike'] == x]           # new df where you are only looking at entries associated with that spiked amount
        
        data['Workflow'] = workflow
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



        if len(current) > 0:                        #If the peptide is found at the spike level
            found_at_spike += 1
            grouped = current.groupby('Rep', as_index= False).agg({'Spike':'first','Modified.Sequence':'first','Precursor.Quantity':'sum' })      #Condense to one line per replicate, sum quantities
            num_reps = len(grouped['Rep'])                                                                                                        #Number of replicates it's found in
            g = mean(grouped['Precursor.Quantity'])
            quant = math.log10(mean(grouped['Precursor.Quantity']))                                                                               #Take the log of the mean quantity
            cv = (grouped['Precursor.Quantity'].std() / mean(grouped['Precursor.Quantity'])) * 100        #Calculate percent CV between quantities found across replicates

            data['Quant_EachRep'] = [x for x in grouped['Precursor.Quantity']]
            data['Number of Reps'] = num_reps
            data['Mean Quant'] = mean(grouped['Precursor.Quantity'])
            data['Log Mean Quant'] = quant
            data['Percent CV'] = cv
            data['Expected Ratio'] = conditions[-4] / x
            # data['Type Phosphopeptide'] = 'Heavy Human'


            #Add reps and CVs information to summarizing table
            summary_heavies[spike_label]['Reps'].append(num_reps)
            summary_heavies[spike_label]['CVs'].append(cv)

        df_list.append(data)

    df = pd.DataFrame(df_list)
    


    if found_at_spike > 0:      #If the peptide is detected at any spike level, this counter should be > 0. Set boolean to true.
        Found = True
        
        info_path = report_directory_path + platform + '_AllSpikes_Info/'
        makePath(info_path)
        x = find.replace(':', '')
        df.to_csv(info_path + x + '_' + workflow + '_output.tsv', sep='\t')
    else:
        Found = False



###HEAVY RATIOS### This is done where each injection a peptide is found in, the quant is ratio'd to the mean quant at the 'set' point
    if Found == True:
        ratios = []
        set = df[df['Spike'] == norm_spike]

        quants = set['Quant_EachRep'].tolist()[0]
        if len(quants) > 0:               #THIS IS WHERE I AM LOSING IDS, some peptides are not detected at 1 fmol, so there's no quant to normalize to

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


            x = find.replace(':', '')
            new_df.to_csv(path + x + '_replicates_output.tsv', sep='\t')
            new_df.to_csv(all_spikes_path + x + '_replicates_output.tsv', sep='\t')

