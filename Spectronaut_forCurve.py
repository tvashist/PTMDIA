import pandas as pd
import re
import math
from statistics import mean,stdev,mode

##Hey!!!!
#Uncomment which instrument you're working on data from
# platform = 'Exploris'
platform = 'Pro'
# platform = 'SCP'

if platform == 'SCP':
    conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0]               #Spiked amounts of peptide used for this instrument platform
    spectronaut = pd.read_csv('S:/Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SCP20220618_114730_PTMDIAProject_SCP_PhosphoBGCurve_Report.tsv', delimiter='\t')



if platform == 'Pro':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv('S:/Tan/PTMDIAProject_PhosphoBGCurve/Outputs/Pro/20220714_100348_PTMDIAProject_Pro_PhosphoBG_Report.tsv', delimiter= '\t')


if platform == 'Exploris':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv('', delimiter = '\t')


###Lights###
summary_lights = {}
lights = pd.read_csv('Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut

for sequence in lights:


    find = sequence     #Spiked in peptide sequence


    single = spectronaut.loc[spectronaut['EG.IntPIMID'] == find]            #Only look at the entries where that specific peptide was found
    single = single[['R.Condition','R.Replicate','EG.IntPIMID','FG.Quantity']]      #Keep only the column indicating condition/spike, peptide sequence, and quantity
    single = single[single['FG.Quantity'].notnull()]                        #Remove columns where quantity is null

    #Add new column with spiked amount
    spiked = []
    for x in single['R.Condition']:
        amount = float(x.replace('fmol',''))
        spiked.append(amount)

    single['Spiked'] = spiked


    df_list = []
    for x in conditions:                                    #For every point on the curve, create a new data frame

        spike_label = str(x) + 'fmol'
        if spike_label not in summary_lights:               #If the peptide isn't already in the dictionary, add
            summary_lights[spike_label] = {}
            summary_lights[spike_label]['Reps'] = []
            summary_lights[spike_label]['CVs'] = []

        data = {}
        current = single[single['Spiked'] == x]             #New df where you are only looking at entries associated with that spiked amount

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

        #If spiked peptide is found at this point on the curve
        if len(current) > 0:
            grouped = current.groupby('R.Replicate', as_index= False).agg({'R.Condition':'first','Spiked':'first','EG.IntPIMID':'first','FG.Quantity':'sum' })      #Condense to one line per replicate, sum quantities
            num_reps = len(grouped['R.Replicate'])              #Number of replicates it's found in
            quant = math.log10(mean(grouped['FG.Quantity']))    #Take the log of the mean quantity
            cv = (grouped['FG.Quantity'].std() / mean(grouped['FG.Quantity'])) * 100        #Calculate percent CV between quantities found across replicates

            #Add relevant data to dictionary
            data['Number of Reps'] = num_reps
            data['Mean Quant'] = mean(grouped['FG.Quantity'])
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
        df.to_csv('PhosphoBG_Curve/' +platform+'_Lights_outputs_NotFound/' + find + '_output.tsv', sep='\t')
        # df.to_csv('Z:/LabMembers/Tan/DIA_QuantitativePTMs/Phospho_TestCurve/Found_Lights_heavies_Curves/' + platform + '_Lights_outputs_NotFound/' + find + '_output.tsv', sep='\t')

    else:                                   #If peptide is found at any point in the curve, send it to the 'Found' folder
        df.to_csv('PhosphoBG_Curve/' +platform+'_Lights_outputs_Found/' + find+'_output.tsv', sep = '\t')
        # df.to_csv('Z:/LabMembers/Tan/DIA_QuantitativePTMs/Phospho_TestCurve/Found_Lights_heavies_Curves/' +platform + '_Lights_outputs_Found/' + find + '_output.tsv', sep='\t')

summary_lights_compiled = {}
for point in summary_lights:
    summary_lights_compiled[point] ={}
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
summary.to_csv('PhosphoBG_Curve/' +platform+'_Lights_outputs_Found/' + platform+'_Lights_Summary.tsv', sep = '\t')






###Heavies###
#All the same comments applied to the lights code above apply to heavies below
summary_heavies = {}

heavies = pd.read_csv('Modified_Heavies.tsv', delimiter= '\t')['Modified']
found_heavies = spectronaut.loc[spectronaut['EG.PTMLocalizationProbabilities'].str.contains('Label')]        #Only look at the entries that contain heavy modification


for sequence in heavies:
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
        data['Spike'] = str(x) +'fmol'
        data['Log Spike'] = math.log10(x)

        if len(current) == 0:
            data['Number of Reps'] = 0
            data['Mean Quant'] = None
            data['Log Mean Quant'] = None
            data['Percent CV'] = None
            data['Expected Ratio'] = conditions[-4] / x


        if len(current) > 0:
            grouped = current.groupby('R.Replicate', as_index= False).agg({'R.Condition':'first','Spiked':'first','EG.IntPIMID':'first','FG.Quantity':'sum' })
            num_reps = len(grouped['R.Replicate'])
            quant = math.log10(mean(grouped['FG.Quantity']))
            cv = (grouped['FG.Quantity'].std() / mean(grouped['FG.Quantity'])) * 100

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


    if (df['Number of Reps'] == 0).all():
        df.to_csv('PhosphoBG_Curve/' + platform+'_Heavies_outputs_NotFound/' + find + '_output.tsv', sep='\t')
        # df.to_csv('Z:/LabMembers/Tan/DIA_QuantitativePTMs/Phospho_TestCurve/Found_Lights_heavies_Curves/' + platform + '_Heavies_outputs_NotFound/' + find + '_output.tsv',sep='\t')

    else:
        df.to_csv('PhosphoBG_Curve/' +platform+'_Heavies_outputs_Found/' + find+'_output.tsv', sep = '\t')
        # df.to_csv('Z:/LabMembers/Tan/DIA_QuantitativePTMs/Phospho_TestCurve/Found_Lights_heavies_Curves/' + platform + '_Heavies_outputs_Found/' + find + '_output.tsv', sep='\t')


summary_heavies_compiled = {}
for point in summary_heavies:
    summary_heavies_compiled[point] ={}
    summary_heavies_compiled[point]['Mode # Reps'] = mode(summary_heavies[point]['Reps'])

    #frequency at which any peptide is seen in 3,2,1 replicates for that given point on the curve
    summary_heavies_compiled[point]['Seen 3x'] = (summary_heavies[point]['Reps']).count(3)
    summary_heavies_compiled[point]['Seen 2x'] = (summary_heavies[point]['Reps']).count(2)
    summary_heavies_compiled[point]['Seen 1x'] = (summary_heavies[point]['Reps']).count(1)

    filtered_cvs = [x for x in summary_heavies[point]['CVs'] if math.isnan(x) == False]          #remove nan values to calculate mean cvs
    summary_heavies_compiled[point]['Mean CV'] = mean(filtered_cvs)

summary = pd.DataFrame.from_dict((summary_heavies_compiled))
summary.to_csv('PhosphoBG_Curve/' +platform+'_Heavies_outputs_Found/' + platform+'_Heavies_Summary.tsv', sep = '\t')









