import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


###USER INPUT###

#Uncomment which instrument you're working on data from
# platform = 'Exploris'
# platform = 'Pro_SmallLibSearch_LocFilter'
# platform = 'SCP'
platform = 'Pro_12fxnOnlySearch_LocFilter'

report_directory_path = 'Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/' + platform + "/"     #Where is your Spectronaut output report?


#Spectronaut output reports found in helium
if platform == 'SCP':
    conditions = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0]               #Spiked amounts of peptide used for this instrument platform
    spectronaut = pd.read_csv(report_directory_path + '20220813_024041_PTMDIAProject_SCP_DIACurveAnalysis_WithSpecLib_Report.tsv', delimiter='\t', low_memory = False)


if platform == 'Pro_SmallLibSearch_LocFilter':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220902_105134_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_SmallLib0.75Loc_Report .tsv', on_bad_lines= 'skip', delimiter= '\t',low_memory = False)


if platform == 'Pro_12fxnOnlySearch_LocFilter':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220906_093527_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_12fxnLib0.75Loc_Report.tsv',delimiter='\t', low_memory=False)
    print("READ")

if platform == 'Pro':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv(report_directory_path + '20220803_134456_PTMDIAProject_DIACurveAnalysis_WithSpecLib_Report_addedFGLabel.tsv', on_bad_lines= 'skip', delimiter= '\t',low_memory = False)


if platform == 'Exploris':
    conditions = [0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
    spectronaut = pd.read_csv('', delimiter = '\t')



#
###Lights###
summary_lights = {}

#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut

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
        data['Spike'] = str(x) + 'fmol'
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

    if (df['Number of Reps'] == 0).all():   #If peptide is not found on any point on the curve, it is not found at all\
        path = report_directory_path + platform + '_Lights_outputs_NotFound/'

        #Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        df.to_csv(path + find + '_output.tsv', sep='\t')




    else:                                   #If peptide is found at any point in the curve, send it to the 'Found' folder
        path = report_directory_path + platform + '_Lights_outputs_Found/'

        # Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        df.to_csv(path + find+'_output.tsv', sep = '\t')


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






###Heavies###
#All the same comments applied to the lights code above apply to heavies below
summary_heavies = {}

heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Heavies.tsv', delimiter= '\t')['Modified_HeaviesAnnotated'][0:234]
found_heavies = spectronaut.loc[spectronaut['FG.LabeledSequence'].str.contains('Label')]        #Only look at the entries that contain heavy modification


new_annotations = []
for seq in found_heavies['FG.LabeledSequence']:      #This column contains sequences with all mod annotations including heavy labeling
    new = seq.replace('Phospho (STY)', '+80').replace('Carbamidomethyl (C)','+57').replace('Label:13C(6)15N(4)','+10').replace('Label:13C(6)15N(2)','+8').replace('[Oxidation (M)]', '').replace('[Label:13C(5)15N(1)]', '').replace('[Label:13C(6)15N(2)]','').replace('[Label:T-5]','').replace('[Label:13C(6)15N(1)]','').replace('[Acetyl (Protein N-term)]','')
    new = new[1:-1]
    new_annotations.append(new)
found_heavies['AnnotatedSequence'] = new_annotations

for sequence in heavies:
    find = sequence     #Spiked in peptide sequence


    single = found_heavies.loc[found_heavies['AnnotatedSequence'] == find]        #Only look at the entries where that specific peptide was found
    single = single[['R.Condition','R.Replicate','FG.LabeledSequence','FG.Quantity', 'AnnotatedSequence']]
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
            grouped = current.groupby('R.Replicate', as_index= False).agg({'R.Condition':'first','Spiked':'first','FG.LabeledSequence':'first','FG.Quantity':'sum' })
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
        else:
            error.append(None)

    df['Percent Error'] = error     #Add new column to dataframe


    if (df['Number of Reps'] == 0).all():
        path = report_directory_path + platform + '_Heavies_outputs_NotFound/'

        #Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        df.to_csv(path + str(find) + '_output.tsv', sep='\t')

    else:
        path = report_directory_path + platform + '_Heavies_outputs_Found/'

        # Create the directory if it doesn't already exist
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

        df.to_csv(path + find + '_output.tsv', sep='\t')

summary_heavies_compiled = {}
for point in summary_heavies:
    summary_heavies_compiled[point] ={}
    if len(summary_heavies[point]['Reps']) != 0:
        summary_heavies_compiled[point]['Mode # Reps'] = mode(summary_heavies[point]['Reps'])

    #frequency at which any peptide is seen in 3,2,1 replicates for that given point on the curve
    summary_heavies_compiled[point]['Seen 3x'] = (summary_heavies[point]['Reps']).count(3)
    summary_heavies_compiled[point]['Seen 2x'] = (summary_heavies[point]['Reps']).count(2)
    summary_heavies_compiled[point]['Seen 1x'] = (summary_heavies[point]['Reps']).count(1)

    filtered_cvs = [x for x in summary_heavies[point]['CVs'] if math.isnan(x) == False]          #remove nan values to calculate mean cvs
    if len(filtered_cvs) != 0:
        summary_heavies_compiled[point]['Mean CV'] = mean(filtered_cvs)

summary = pd.DataFrame.from_dict((summary_heavies_compiled))

heavies_exist_path = report_directory_path + platform + '_Heavies_outputs_Found/'
heaviesExist = os.path.exists(heavies_exist_path)
if heaviesExist:
    summary.to_csv(report_directory_path + platform + '_Heavies_outputs_Found/' + platform + '_Heavies_Summary.tsv', sep = '\t')









