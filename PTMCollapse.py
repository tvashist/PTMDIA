import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


wd = 'Y:\\LabMembers\\Tan\\CompGroup\\PTMSiteCollapse\\'
# file = pd.read_csv(wd + '20221122_163329_PTMDIAProject_ExplorisFAIMS_3SSLib_DIACurveAnalysis_Report.tsv', delimiter = '\t', low_memory= False)



report = pd.read_csv(wd + 'TwoFiles.tsv', delimiter ='\t', low_memory= False)

# Only keep sequences that contain the target modification ([+80])
target = report.loc[report['EG.IntPIMID'].str.contains('\[\+80\]')]

# Create new column that removes other modifications besides mod of interest, and new column with modified site position information
def find_position(target_sequence):     #Finds position of ACTUALLY MODIFIED sty's in peptide sequence, and also returns list of annotations like in SM (i.e. S18, T24)
    locations = []
    annotations = []

    pos = target_sequence.find('[+80]')     #the first phosphosite on the sequence

    while pos != -1:
        aa = target_sequence[pos-1]
        annotations.append(aa + str(pos))

        locations.append(pos)
        target_sequence = target_sequence[:pos] + target_sequence[pos+5:]
        pos = target_sequence.find('[+80]')


    return(annotations,locations) #returns residues that are actually modified

# print(find_position('SSQTGTEPGAAHTSS[+80]PQPSTSR'))

def sty_ProteinLocations(PTM_ProteinLocations): #Because of this function, miscleavages will still have the same collapse key, sequences with different non-STY mods will have the same key
    sty_locations = []

    isoforms = PTM_ProteinLocations.split(';')
    for iso in isoforms:
        sites = iso.replace('(','').replace(')','').split(',')

        sty_only = []
        for s in sites:
            if s.startswith(('S','T','Y')):
                sty_only.append(s)

        sty_locations.append(sty_only)

    return(sty_locations)

# print(sty_ProteinLocations('(M690,S691,S693);(M662,S663,S665);(M97,S98,S100);(M711,S712,S714)'))


target_mods = []
actually_modified = []
PTM_PeptideLocations = []
collapse_keys = []

#Append columns to existing data frame, including collapse key and other annotations
for index, row in target.iterrows():
    all_mods = row['EG.IntPIMID']       #This value annotates all mods

    target_mods_seq = all_mods.replace('[+16]','').replace('[+57]','').replace('[+42]','').replace('_','')  #Remove M-oxidation, carbamidomethylation, and N-term acetylation
    target_mods.append(target_mods_seq)   #Add altered sequence to new column

    position_info = find_position(target_mods_seq)
    PTM_PeptideLocations.append(position_info[0])           #[S18, T24]
    actually_modified.append(position_info[1])              #[18, 24]

    protein_ids = row['PG.ProteinAccessions']
    PTM_ProteinLocations = row['EG.ProteinPTMLocations']                                #Existing column also contains PTM locations M-Ox, Cam, etc
    STY_ProteinLocations = str(sty_ProteinLocations(PTM_ProteinLocations))              #Create column that only contains PTM locations of STY's, this way even sequences with other modifications will still collapse if they have the same phosphorylation

    collapse_keys.append((STY_ProteinLocations, protein_ids))

target['Only_Target_Mods'] = target_mods
target['PTM_PeptidePositions'] = actually_modified
target['PTM_PeptidePositions_AA'] = PTM_PeptideLocations
target['Collapse'] = collapse_keys

#Create new column that sums the quantity across all grouped rows; that way when selecting the representative column, the sum will still be in one column
target['FG.SummedQuants'] = target.groupby(['Collapse'])['FG.Quantity'].transform('sum')
# target.to_csv(wd + 'NewSumColumn.tsv', sep= '\t', index = False)

#Only keep row from a group that contains the highest PTM Assay Probability Score (most confident localization); There will be repeats for rows with same PTM Assay Probability Score
group_MaxPTMAssayScore = target.groupby(['Collapse'])['EG.PTMAssayProbability'].transform(max) == target['EG.PTMAssayProbability']
collapsed_MaxPTMAssayScore = target[group_MaxPTMAssayScore]

#Make collapse key first column, will remove this helper column for final output
first_column = collapsed_MaxPTMAssayScore.pop('Collapse')
collapsed_MaxPTMAssayScore.insert(0, 'Collapse', first_column)

#At this stage, all repeated lines have the same PTM assay probability score

#Only keep row from group with highest EG.Cscore
group_MaxCScore = collapsed_MaxPTMAssayScore.groupby(['Collapse'])['EG.Cscore'].transform(max) == collapsed_MaxPTMAssayScore['EG.Cscore']
collapsed_MaxCScore = collapsed_MaxPTMAssayScore[group_MaxCScore]

#At this stage, there should not be replicates



collapsed_MaxCScore.to_csv(wd + 'Collapsed_CScore.tsv', sep= '\t', index = False)






















