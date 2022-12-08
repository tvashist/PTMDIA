import pandas as pd
pd.options.mode.chained_assignment = None
import re
import math
from statistics import mean,stdev,mode
import os
import numpy as np

wd = 'Z:\\Alexi\\TimsTOF_Jurkat_iRTs_phospho\\Results\\PTM_results\\'
file_name = '20221205_153854_Phospho_Jurkat_with_Pro12FxLibv1_Report_PEPPositions.tsv'



report = pd.read_csv(wd + file_name, delimiter ='\t', low_memory= False)
# short = report[:100]
# short.to_csv(wd + 'Toy.tsv', sep = '\t', index = False)

print('File Read')


samples = (report['R.FileName']).unique()       #List of unique samples that IDs may come from


#Function works without bugs
def find_heavies(reg, heavy):               #Inputs are non-heavy modified sequence and heavy modified version of sequence
    has_heavies = True

    reg_split = re.split('\[|\]',reg)
    heavy_split = re.split('\[|\]',heavy)

    reg_annotations = set([x for x in reg_split if x.startswith('+')])
    heavy_annotations = set([x for x in heavy_split if x.startswith('+')])


    diff = heavy_annotations - reg_annotations              #This will give list of unique heavy mods that are present in this entry

    collect = []
    for element in diff:
        expression = '\\[\\' + element + '\\]'              #Generates regex for each heavy mod
        collect.append(expression)

    if len(diff) == 0:                                      #If there are no heavy mods, return None
        has_heavies = False

    #If heavies are present
    if len(diff) == 1:              #Only one heavy mod type
        prep_re = collect[0]

    if len(diff) >1:                #More than one heavy mod type
        prep_re = '|'.join(collect)

    # print(type(prep_re))

    if has_heavies == True:         #Returns regular expression if heavies are present
        return(prep_re)

    if has_heavies == False:
        return(None)


# print(find_heavies('EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK'))






def remove_experimental_mods(sequence):
    new = sequence.replace('[+16]', '').replace('[+57]', '').replace('[+42]', '').replace('_', '')      #Remove carbamidomethylation, N-term acetylation, M-oxidation
    return(new)

# Only keep sequences that contain the target modification ([+80])
target = report.loc[report['FG.IntMID'].str.contains('\[\+80\]')]

def find_position(non_heavy_mods,target_sequence, PEP_position):     #Finds position of ACTUALLY MODIFIED sty's in peptide sequence, and also returns list of annotations like in SM (i.e. S18, T24); ALSO FIND POSITION OF HEAVY MODS
    #After finding the peptide positions of STY, remove [+80] annotations so you can focus on just the heavy annotations remaining
    heavies = target_sequence.replace('[+80]', '')

    has_heavies = False
    heavy_ProteinPositions = []
    heavy_annotations = []
    heavy_locations = []

    chars = []
    pos = []

    exp = find_heavies(non_heavy_mods,target_sequence)
    if exp != None:
        has_heavies = True
        finder = re.finditer(exp, heavies)
        for match in finder:
            chars.append(match.group())
            pos.append(match.start())



    while len(pos) > 0:
        h_pos = min(pos)
        aa = heavies[h_pos -1]

        index_min = np.argmin(pos)
        h_char = chars[index_min]


        heavy_annotations.append(aa + h_char)
        heavy_locations.append(h_pos)

        skip = len(h_char)
        heavies = heavies[:h_pos] + heavies[h_pos + skip:]

        chars = []
        pos = []

        finder = re.finditer(exp, heavies)
        for match in finder:
            chars.append(match.group())
            pos.append(match.start())

    for x in PEP_position.split(';'):
        annotation = []
        for y in range(0, len(heavy_locations)):
            prot_pos = int(x) + heavy_locations[y] - 1
            annotation.append(heavy_annotations[y] + str(prot_pos))
        heavy_ProteinPositions.append(annotation)

    if has_heavies == True:
        return(heavy_ProteinPositions) #returns residues that are actually modified
    else:
        return(None)
# print(find_position('EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','EGPTDHLESACPL[+8]NLPL[+10]QNNHT[+6]AADMYLS[+80]PVRS[+80]PK','580;570'))
# print(find_position('EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','580;570'))

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
sty_PeptideLocations = []
heavy_ProteinLocations = []
collapse_keys = []

for index, row in target.iterrows():
    all_mods = row['FG.IntMID']       #This value annotates all mods
    non_heavy_mods = row['EG.IntPIMID']


    PEP_pos = row['PEP.PeptidePosition']

    target_mods_seq = remove_experimental_mods(all_mods)  #Remove M-oxidation, carbamidomethylation, and N-term acetylation, INCLUDE HEAVIES
    only_targets_seq = remove_experimental_mods(non_heavy_mods)

    target_mods.append(only_targets_seq)   #Add altered sequence to new column

    position_info = find_position(non_heavy_mods, target_mods_seq, PEP_pos)
    # sty_PeptideLocations.append(position_info[0])           #[S18, T24]
    heavy_ProteinLocations.append(position_info)         #[L[+10]564]


    protein_ids = row['PG.ProteinAccessions']
    PTM_ProteinLocations = row['EG.ProteinPTMLocations']                                #Existing column also contains PTM locations M-Ox, Cam, etc
    STY_ProteinLocations = str(sty_ProteinLocations(PTM_ProteinLocations))              #Create column that only contains PTM locations of STY's, this way even sequences with other modifications will still collapse if they have the same phosphorylation

    k = (protein_ids,STY_ProteinLocations,str(position_info))
    collapse_keys.append(k)



target['Only_Target_Mods'] = target_mods
# target['PTM_PeptidePositions'] = sty_PeptideLocations
target['HeavyMod_ProteinPositions'] = heavy_ProteinLocations
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


keep = ['Collapse','R.Condition','R.FileName', 'R.Replicate', 'PG.Genes','PG.Organisms', 'PG.ProteinAccessions', 'PG.IBAQ', 'PG.Quantity', 'PEP.PeptidePosition', 'PEP.Quantity', 'EG.IntPIMID', 'EG.ProteinPTMLocations','EG.PTMPositions [Phospho (STY)]',
        'EG.PTMAssayProbability','EG.PTMLocalizationProbabilities','EG.Cscore','EG.IntCorrScore','FG.Charge','FG.IntMID','FG.CalibratedMassAccuracy (PPM)', 'Only_Target_Mods', 'HeavyMod_ProteinPositions']

for s in samples:
    keep.append(s)


collapsed_report = collapsed_MaxCScore[keep]

def summarize_any(full, collapsed_report, sample):
    all_sequences = full['FG.IntMID']  # All precursors
    num_modified_sequences = len(
        all_sequences.unique())  # All precursors, including experimentally-introduced mods and heavy mods

    phos_sequences = [x for x in all_sequences if '[+80]' in x]  # Only phospho-containing precursors
    enrichment = (len(set(phos_sequences)) / num_modified_sequences) * 100  # Phospho-containing precursors / all precursors

    # Number of unique STY modified sequences (Heavy mods allowed, no experimentally-introduced mods)
    phospho_no_experimental = list(map(remove_experimental_mods, phos_sequences))  # Remove experimentally-introduced modifications
    unique_phosphopeptides = len(set(phospho_no_experimental))  # This will represent the number of phosphopeptides

    if sample != 'Combined':
        select = collapsed_report[sample].values.tolist()
        res = [i for i in select if i is not None]
        phosphosites = len(res)

    if sample == 'Combined':
        phosphosites = len(collapsed_report)

    return({'Run': sample, '%_Precursors_with_sty': enrichment, 'num_modified_sequences': num_modified_sequences ,
                            'num_phosphopeptides': unique_phosphopeptides, 'num_phosphosites': phosphosites})

def summarize(report, collapsed_report, samples):
    print('Summarizing data...')
    data = pd.DataFrame(columns= ['Run','%_Precursors_with_sty', 'num_modified_sequences', 'num_phosphopeptides', 'num_phosphosites'])       #Initialize datarframe with columns

    for s in samples:
        one = report[report['R.FileName'] == s]
        ret = summarize_any(one,collapsed_report, s)
        data = data.append(ret, ignore_index= True)

    total = summarize_any(report, collapsed_report, 'Combined')
    data = data.append(total, ignore_index= True)


    data.to_csv(wd + 'VMSummary_02.tsv', sep = '\t', index = False)



summarize(report, collapsed_report, samples)
collapsed_report.to_csv(wd + 'VMReport_CollapseKeys_AddedFunc.tsv', sep= '\t', index = False)



















 ###########################################
    # #can't find a work around to look for first instance of either-or
    #
    # finder = min(("[+8]", heavies.find('[+8]')), ("[+10]", heavies.find('[+10]')), key=lambda x: x[1])
    # if finder[1] == -1:
    #     finder = max(("[+8]", heavies.find('[+8]')), ("[+10]", heavies.find('[+10]')), key=lambda x: x[1])
    #
    # if finder[1] == -1:         #If there are no heavy amino acids in the sequence
    #     h_pos = -1
    #     has_heavies = False
    # else:
    #     has_heavies = True
    #     h_pos = finder[1]
    #     h_char = finder[0]
    #
    # while h_pos != -1:
    #     aa = heavies[h_pos-1]                       #The heavy labeled amino acid
    #
    #     heavy_annotations.append(aa+h_char)         #Collect annotation of the amino acid and the heavy label (+8/+10)
    #     heavy_locations.append(h_pos)               #Position of heavy-labeled amino acid within the (stripped) peptide sequence
    #     skip = len(h_char)                          #How long is the heavy annotation in the sequence ([+8] has 4 chars rather than 5)
    #
    #     heavies = heavies[:h_pos] + heavies[h_pos+skip:]
    #     # print(heavies)
    #
    #     finder = min(("[+8]", heavies.find('[+8]')), ("[+10]", heavies.find('[+10]')), key=lambda x: x[1])
    #     if finder[1] == -1:
    #         finder = max(("[+8]", heavies.find('[+8]')), ("[+10]", heavies.find('[+10]')), key=lambda x: x[1])
    #
    #
    #     h_pos = finder[1]
    #     h_char = finder[0]
############################################################