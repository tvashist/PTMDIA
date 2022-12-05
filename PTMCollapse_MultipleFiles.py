import pandas as pd
pd.options.mode.chained_assignment = None
import re
import math
from statistics import mean,stdev,mode
import os

wd = 'Z:\\Alexi\\TimsTOF_Jurkat_iRTs_phospho\\Results\\PTM_results\\'
file_name = '20221205_153854_Phospho_Jurkat_with_Pro12FxLibv1_Report_PEPPositions.tsv'


report = pd.read_csv(wd + file_name, delimiter ='\t', low_memory= False)
print('File Read')


samples = (report['R.FileName']).unique()       #List of unique samples that IDs may come from



# Only keep sequences that contain the target modification ([+80])
target = report.loc[report['EG.IntPIMID'].str.contains('\[\+80\]')]

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


target_mods = []
actually_modified = []
PTM_PeptideLocations = []
collapse_keys = []

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



#In this section, each group under a collapse key is further split into members associated with each individual sample, then quant is separately summed per run within the gorup and added to a dictionary
#Each collapse key is a key in a dictionary where the values

groups = target.groupby('Collapse')

meta = {}
for name, group in groups:
    meta[name] = None

    file_quant = {}

    for s in samples:           #Sum quant by file
        file_quant[s] = None

        sub = group[group['R.FileName'] == s]
        if len(sub) > 0:
            quant = sub['FG.Quantity'].sum()
        else:
            quant = None

        file_quant[s] = quant

    meta[name] = file_quant


#Only keep row from a group that contains the highest PTM Assay Probability Score (most confident localization); There will be repeats for rows with same PTM Assay Probability Score
group_MaxPTMAssayScore = target.groupby(['Collapse'])['EG.PTMAssayProbability'].transform(max) == target['EG.PTMAssayProbability']
collapsed_MaxPTMAssayScore = target[group_MaxPTMAssayScore]

#Make collapse key first column, will remove this helper column for final output
first_column = collapsed_MaxPTMAssayScore.pop('Collapse')
collapsed_MaxPTMAssayScore.insert(0, 'Collapse', first_column)

#Only keep row from group with highest EG.Cscore
group_MaxCScore = collapsed_MaxPTMAssayScore.groupby(['Collapse'])['EG.Cscore'].transform(max) == collapsed_MaxPTMAssayScore['EG.Cscore']
collapsed_MaxCScore = collapsed_MaxPTMAssayScore[group_MaxCScore]



for f in samples:
    collapsed_MaxCScore[f] = [None for x in range(0,len(collapsed_MaxCScore))]      #Create new exmpty column for quant by sample


for index, row in collapsed_MaxCScore.iterrows():          #Iterate through each entry
    key = row['Collapse']


    for f in samples:                                       #Add quant for this row based on sample it's coming from
        q = meta[key][f]
        collapsed_MaxCScore.at[index ,f] = q


keep = ['R.Condition','R.FileName', 'R.Replicate', 'PG.Genes','PG.Organisms', 'PG.ProteinAccessions', 'PG.IBAQ', 'PG.Quantity', 'PEP.PeptidePosition', 'PEP.Quantity', 'EG.IntPIMID', 'EG.ProteinPTMLocations','EG.PTMPositions [Phospho (STY)]',
        'EG.PTMAssayProbability','EG.PTMLocalizationProbabilities','EG.Cscore','EG.IntCorrScore','FG.Charge','FG.IntMID','FG.CalibratedMassAccuracy (PPM)', 'Only_Target_Mods', 'PTM_PeptidePositions_AA']

for s in samples:
    keep.append(s)

for s in samples:
    select = collapsed_MaxCScore[s].values.tolist()
    res = [i for i in select if i is not None]
    print(s)
    print(len(res))




report = collapsed_MaxCScore[keep]
report.to_csv(wd + 'VMReport.tsv', sep= '\t', index = False)

