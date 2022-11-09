from unittest import skip

import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


#Uncomment which instrument you're working on data from
# platform = 'Exploris'
platform = 'Exploris_FAIMS'
# platform = 'TimsTOF_Pro'
# platform = 'TimsTOF_SCP'

#Library
# library = 'Library_3SS_Spiked'
# library = 'Library_12fxn_NotSpiked'
library = 'Library_Combined'


report_directory_path = 'Z:/Helium_Tan/FINAL_PTMDIA/' + platform + '/DIANN/' + library + '/SearchOutputs/dia_nn/out/'
diann = pd.read_csv(report_directory_path+ 'report.tsv', delimiter='\t', low_memory=False)
print("Report read")

if platform == 'TimsTOF_SCP':
    conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0]               #Spiked amounts of peptide used for this instrument platform

else:
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]


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
        rep = int(parts[-1])

    run.append(rep)

diann['Rep'] = run
diann['Spike'] = spike


# diann = diann[diann['Rep'] != 0]            #Only look at rows without spike



# ###Lights###
summary_lights = {}

#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['DIANN_Mods']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut



for sequence in lights:


    find = str(sequence)     #Spiked in peptide sequence


    single = diann.loc[diann['Modified.Sequence'] == find]            #Only look at the entries where that specific peptide was found
    single = single[['Rep','Spike','Modified.Sequence','Precursor.Quantity']]      #Keep only the column indicating condition/spike, peptide sequence, and quantity
    single = single[single['Precursor.Quantity'].notnull()]                        #Remove columns where quantity is null

    df_list = []
    for x in conditions:                                    #For every point on the curve, create a new data frame

        spike_label = str(x) + 'fmol'
        if spike_label not in summary_lights:               #If the peptide isn't already in the dictionary, add
            summary_lights[spike_label] = {}
            summary_lights[spike_label]['Reps'] = []
            summary_lights[spike_label]['CVs'] = []

        data = {}
        current = single[single['Spike'] == x]             #New df where you are only looking at entries associated with that spiked amount
        data['Type Phosphopeptide'] = 'Light Non-Human'     #Annotate whether we spiked in non-human light or human heavy
        data['Peptide'] = find
        data['Spike'] = spike_label
        data['Log Spike'] = math.log10(x)                   #Plot log spiked amount vs log quantity

        #If spiked peptide is not found at this point on the curve
        if len(current) == 0:
            data['Number of Reps'] = 0
            data['Mean Quant'] = None
            data['Log Mean Quant'] = None
            data['Percent CV'] = None
            data['Expected Ratio'] = conditions[-4] / x
#
        #If spiked peptide is found at this point on the curve
        if len(current) > 0:
            grouped = current.groupby('Rep', as_index= False).agg({'Spike':'first','Modified.Sequence':'first','Precursor.Quantity':'sum' })      #Condense to one line per replicate, sum quantities

            num_reps = len(grouped['Rep'])              #Number of replicates it's found in
            quant = math.log10(mean(grouped['Precursor.Quantity']))    #Take the log of the mean quantity
            cv = (grouped['Precursor.Quantity'].std() / mean(grouped['Precursor.Quantity'])) * 100        #Calculate percent CV between quantities found across replicates

            #Add relevant data to dictionary
            data['Number of Reps'] = num_reps
            data['Mean Quant'] = mean(grouped['Precursor.Quantity'])
            data['Log Mean Quant'] = quant
            data['Percent CV'] = cv
            data['Expected Ratio'] = conditions[-4] / x

            #Add reps and CVs information to summarizing table
            summary_lights[spike_label]['Reps'].append(num_reps)
            summary_lights[spike_label]['CVs'].append(cv)


        df_list.append(data)        #Collect all dictionaries into a list, so by the end it will be a list of dicts for each point on the curve


    df = pd.DataFrame(df_list)      #Combine all dicts into one data frame and export separately for each spiked peptide


    #collect ratios of quantities relative to a set point. 1.0fmol for Pro/Exploris and 0.2fmol for SCP
    ratios = []
    set = list(df['Mean Quant'])[-4]
    if set != None:
        for index, row in df.iterrows():
            current = row['Mean Quant']
            if current != None:
                rat = set / current
                ratios.append(rat)
    else:
        for index, row in df.iterrows():
            ratios.append(None)

    df['Actual Ratios'] = ratios

    #Calculates percent error between expected and actual ratios of quantities to set point
    error = []
    for index, row in df.iterrows():
        expected = row['Expected Ratio']
        actual = row['Actual Ratios']

        if expected != None and actual != None:
            perc_error = (abs(actual-expected) / expected) *100
            error.append(perc_error)
            # print(perc_error)
        else:
            error.append(None)

    df['Percent Error'] = error     #Add new column to dataframe

    if (df['Number of Reps'] == 0).all():   #If peptide is not found on any point on the curve, it is not found at all
        path = report_directory_path + platform + '_Lights_outputs_NotFound/'

        #Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        # x = re.sub(r'\([^)]*\)', '', find)
        x = find.replace(':', '')
        df.to_csv(path + x + '_output.tsv', sep='\t')




    else:                                   #If peptide is found at any point in the curve, send it to the 'Found' folder
        path = report_directory_path + platform + '_Lights_outputs_Found/'

        # Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        # x = re.sub(r'\([^)]*\)', '', find)
        x = find.replace(':','')
        df.to_csv(path + x +'_output.tsv', sep = '\t')


summary_lights_compiled = {}

for point in summary_lights:
    summary_lights_compiled[point] = {}
    if len(summary_lights[point]['Reps']) != 0:
        summary_lights_compiled[point]['Mode # Reps'] = mode(summary_lights[point]['Reps'])

    #frequency at which any peptide is seen in 3,2,1 replicates for that given point on the curve
    summary_lights_compiled[point]['Seen 3x'] = (summary_lights[point]['Reps']).count(3)
    summary_lights_compiled[point]['Seen 2x'] = (summary_lights[point]['Reps']).count(2)
    summary_lights_compiled[point]['Seen 1x'] = (summary_lights[point]['Reps']).count(1)


    filtered_cvs = [x for x in summary_lights[point]['CVs'] if math.isnan(x) == False]          #remove nan values to calculate mean cvs
    if len(filtered_cvs) > 0:
        summary_lights_compiled[point]['Mean CV'] = mean(filtered_cvs)
    else:
        summary_lights_compiled[point]['Mean CV'] = None


summary = pd.DataFrame.from_dict((summary_lights_compiled))

lights_exist_path = report_directory_path + platform + '_Lights_outputs_Found/'
lightsExist = os.path.exists(lights_exist_path)
if lightsExist:
    summary.to_csv(report_directory_path + platform + '_Lights_outputs_Found/' + platform + '_Lights_Summary.tsv', sep = '\t')


###HEAVIES###

heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['DIANN_Mods']

summary_heavies = {}

counts = []
for sequence in heavies:


    find = str(sequence)     #Spiked in peptide sequence


    single = diann.loc[diann['Modified.Sequence'] == find]            #Only look at the entries where that specific peptide was found
    single = single[['Rep','Spike','Modified.Sequence','Precursor.Quantity']]      #Keep only the column indicating condition/spike, peptide sequence, and quantity
    single = single[single['Precursor.Quantity'].notnull()]                        #Remove columns where quantity is null


    df_list = []
    for x in conditions:                                    #For every point on the curve, create a new data frame

        spike_label = str(x) + 'fmol'
        if spike_label not in summary_heavies:               #If the peptide isn't already in the dictionary, add
            summary_heavies[spike_label] = {}
            summary_heavies[spike_label]['Reps'] = []
            summary_heavies[spike_label]['CVs'] = []

        data = {}
        current = single[single['Spike'] == x]             #New df where you are only looking at entries associated with that spiked amount


        data['Type Phosphopeptide'] = 'Light Non-Human'     #Annotate whether we spiked in non-human light or human heavy
        data['Peptide'] = find
        data['Spike'] = spike_label
        data['Log Spike'] = math.log10(x)                   #Plot log spiked amount vs log quantity

#         #If spiked peptide is not found at this point on the curve
        if len(current) == 0:
            data['Number of Reps'] = 0
            data['Mean Quant'] = None
            data['Log Mean Quant'] = None
            data['Percent CV'] = None
            data['Expected Ratio'] = conditions[-4] / x

        #If spiked peptide is found at this point on the curve
        if len(current) > 0:
            grouped = current.groupby('Rep', as_index= False).agg({'Spike':'first','Modified.Sequence':'first','Precursor.Quantity':'sum' })      #Condense to one line per replicate, sum quantities
            num_reps = len(grouped['Rep'])              #Number of replicates it's found in
            quant = math.log10(mean(grouped['Precursor.Quantity']))    #Take the log of the mean quantity
            cv = (grouped['Precursor.Quantity'].std() / mean(grouped['Precursor.Quantity'])) * 100        #Calculate percent CV between quantities found across replicates
            # counts +=1

            #Add relevant data to dictionary
            data['Number of Reps'] = num_reps
            data['Mean Quant'] = mean(grouped['Precursor.Quantity'])
            data['Log Mean Quant'] = quant
            data['Percent CV'] = cv
            data['Expected Ratio'] = conditions[-4] / x

            #Add reps and CVs information to summarizing table
            summary_heavies[spike_label]['Reps'].append(num_reps)
            summary_heavies[spike_label]['CVs'].append(cv)


        df_list.append(data)        #Collect all dictionaries into a list, so by the end it will be a list of dicts for each point on the curve


    df = pd.DataFrame(df_list)      #Combine all dicts into one data frame and export separately for each spiked peptide


    #collect ratios of quantities relative to a set point. 1.0fmol for Pro/Exploris and 0.2fmol for SCP
    ratios = []
    set = list(df['Mean Quant'])[-4]
    if set != None:
        for index, row in df.iterrows():
            current = row['Mean Quant']
            if current != None:
                rat = set / current
                ratios.append(rat)
    else:
        for index, row in df.iterrows():
            ratios.append(None)

    df['Actual Ratios'] = ratios

    #Calculates percent error between expected and actual ratios of quantities to set point
    error = []
    for index, row in df.iterrows():
        expected = row['Expected Ratio']
        actual = row['Actual Ratios']

        if expected != None and actual != None:
            perc_error = (abs(actual-expected) / expected) *100
            error.append(perc_error)
            # print(perc_error)
        else:
            error.append(None)

    df['Percent Error'] = error     #Add new column to dataframe

    if (df['Number of Reps'] == 0).all():   #If peptide is not found on any point on the curve, it is not found at all
        path = report_directory_path + platform + '_Heavies_outputs_NotFound/'

        #Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        # x = re.sub(r'\([^)]*\)', '', find)
        x = find.replace(':','')
        df.to_csv(path + x + '_output.tsv', sep='\t')





    else:                                   #If peptide is found at any point in the curve, send it to the 'Found' folder


        path = report_directory_path + platform + '_Heavies_outputs_Found/'

        # Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        # x = re.sub(r'\([^)]*\)', '(Mod)', find)
        # counts.append(x)
        x = find.replace(':', '')
        df.to_csv(path + x +'_output.tsv', sep = '\t')


summary_heavies_compiled = {}

for point in summary_heavies:
    summary_heavies_compiled[point] = {}
    if len(summary_heavies[point]['Reps']) != 0:
        summary_heavies_compiled[point]['Mode # Reps'] = mode(summary_heavies[point]['Reps'])

    #frequency at which any peptide is seen in 3,2,1 replicates for that given point on the curve
    summary_heavies_compiled[point]['Seen 3x'] = (summary_heavies[point]['Reps']).count(3)
    summary_heavies_compiled[point]['Seen 2x'] = (summary_heavies[point]['Reps']).count(2)
    summary_heavies_compiled[point]['Seen 1x'] = (summary_heavies[point]['Reps']).count(1)


    filtered_cvs = [x for x in summary_heavies[point]['CVs'] if math.isnan(x) == False]          #remove nan values to calculate mean cvs
    if len(filtered_cvs) > 0:
        summary_heavies_compiled[point]['Mean CV'] = mean(filtered_cvs)
    else:
        summary_heavies_compiled[point]['Mean CV'] = None


summary = pd.DataFrame.from_dict((summary_heavies_compiled))

heavies_exist_path = report_directory_path + platform + '_Heavies_outputs_Found/'
heaviesExist = os.path.exists(heavies_exist_path)
if heaviesExist:
    summary.to_csv(report_directory_path + platform + '_Heavies_outputs_Found/' + platform + '_Heavies_Summary.tsv', sep = '\t')