import pandas as pd
pd.options.mode.chained_assignment = None
import re
import math
from statistics import mean,stdev,mode
import os
import numpy as np

def prepare(wd, file_name):
    report = pd.read_csv(wd + file_name, delimiter ='\t', low_memory= False)        #Read in Spectronaut "Normal Report"
    print('File Read')

    samples = (report['R.FileName']).unique()       #List of unique samples/ LC-MS runs that IDs may come from
    return(report, samples)

def find_heavies(reg, heavy):
    """
    Finds the heavy annotations a sequence contains. The 'EG.IntPIMID' column does not contain heavy annotations such as 'K[+8]' or 'R[+10]',
    but the 'FG.IntPIMID' column does. The only differences between the sequences from each column is the presence of heavy annotations.

    Args:
        reg: row['EG.IntPIMID'], sequence without heavy modification annotations
        heavy: row['FG.IntMID'], sequence with heavy modification annotation if the peptide contains heavy amino acids

    Returns: A regex expression that can be used to search for extracted heavy annotations if present. None is returned if peptide does not contain
    heavy-labeled amino acid.
    """

    has_heavies = True                                      #Boolean capturing whether peptide has heavy-labeled amino acids or not

    #Split sequences by '[' and ']'
    reg_split = re.split('\[|\]',reg)
    heavy_split = re.split('\[|\]',heavy)

    #Extract annotations unique to 'heavy' arg. Ex: '+10' or '+8'
    reg_annotations = set([x for x in reg_split if x.startswith('+')])
    heavy_annotations = set([x for x in heavy_split if x.startswith('+')])


    diff = heavy_annotations - reg_annotations              #This will give list of unique heavy mods that are present in this entry

    collect = []
    for element in diff:
        expression = '\\[\\' + element + '\\]'              #Generates regex for each heavy mod
        collect.append(expression)

    if len(diff) == 0:                                      #If there are no heavy mods, set boolean to False
        has_heavies = False

    #If heavy mods are present
    if len(diff) == 1:              #Only one heavy mod type
        prep_re = collect[0]

    if len(diff) >1:                #More than one heavy mod type, compile regex to search for all
        prep_re = '|'.join(collect)

    # Returns regular expression if heavy annotations are present, return None if not
    if has_heavies == True:
        return(prep_re)

    if has_heavies == False:
        return(None)


# print(find_heavies('EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','EGPTDHLESACPL[+8]NLPLQNNHTAADMYLS[+80]PVR[+10]S[+80]PK'))

def remove_experimental_mods(sequence):
    """
    Removes annotations relating to N-term acetylation, oxidation of methionine, and carbamidomethylation. Also removes '_' flanking sequences to ease indexing.

    Args:
        sequence: Sequence containing experimentally-introduced mod annotations listed above.

    Returns: Peptide sequence without experimentally-introduced mod annotations listed above.
    """
    new = sequence.replace('[+16]', '').replace('[+57]', '').replace('[+42]', '').replace('_', '')      #Remove carbamidomethylation, N-term acetylation, M-oxidation
    return(new)



def find_position(non_heavy_mods,target_sequence, PEP_position):
    """
    Finds positions of heavy-labeled amino acids in every protein ID a sequence is mapped to.

    Args:
        non_heavy_mods: Sequence without heavy annotations, stripped of experimentally-induced mod annotations
        target_sequence: Sequence with heavy annotations, stripped of experimentally-introduced mod annotations
        PEP_position: Position of peptide sequence within protein

    Returns: List of heavy modifications present in a sequence, with amino acid and protein position, for each protein ID the sequence maps to.
    Ex: '[K[+8]562, L[+8]570],[K[+8]570, L[+8]578]]'

    """

    heavies = target_sequence.replace('[+80]', '')  #Remove STY mods so only heavy mod annotations are left

    has_heavies = False                             #Boolean tracks if heavy mods are present

    heavy_ProteinPositions = []
    heavy_annotations = []
    heavy_locations = []

    chars = []
    pos = []

    exp = find_heavies(non_heavy_mods,target_sequence)  #Returns regex expression for the heavies that should be searched for, if present
    if exp != None:
        has_heavies = True
        finder = re.finditer(exp, heavies)              #Find all instances of heavy annotations
        for match in finder:
            chars.append(match.group())                 #Extracts heavy annotation string (Ex: '[+8]', appends to list
            pos.append(match.start())                   #Extracts starting index of match, appends to list



    while len(pos) > 0:             #While there are still mods that haven't been accounted for yet in the sequence
        h_pos = min(pos)            #Only look at first occurrence of a heavy, because this will be the only correct position within the sequence.
        aa = heavies[h_pos -1]      #Which amino acid is heavy-labeled?

        index_min = np.argmin(pos)  #Index of min value in sequence
        h_char = chars[index_min]   #Extracts annotation


        heavy_annotations.append(aa + h_char)
        heavy_locations.append(h_pos)

        skip = len(h_char)          #Length of annotation in sequence
        heavies = heavies[:h_pos] + heavies[h_pos + skip:]  #Remove this annotation as it's now accounted for, generate new sequence with remaining annotations.

        chars = []
        pos = []

        finder = re.finditer(exp, heavies)
        for match in finder:
            chars.append(match.group())
            pos.append(match.start())

    for x in PEP_position.split(';'):               #PEP_positions relating to different protein IDs are separated by ';'
        annotation = []
        for y in range(0, len(heavy_locations)):    #Collect heavy-labeled annotations and positions for each protein ID and respective peptide position
            prot_pos = int(x) + heavy_locations[y] - 1
            annotation.append(heavy_annotations[y] + str(prot_pos))
        heavy_ProteinPositions.append(annotation)

    if has_heavies == True:
        return(heavy_ProteinPositions) #returns residues that are actually modified
    else:
        return(None)
# print(find_position('EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','EGPTDHLESACPL[+8]NLPL[+10]QNNHT[+6]AADMYLS[+80]PVRS[+80]PK','580;570'))
# print(find_position('EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','EGPTDHLESACPLNLPLQNNHTAADMYLS[+80]PVRS[+80]PK','580;570'))

def sty_ProteinLocations(PTM_ProteinLocations):
    """
    Generates part of collapse key without experimentally-introduced mods. Miscleavages, different charge states, sequences with different
    experimentally-introduced mod sites will still be collapsed together.

    Args:
        PTM_ProteinLocations: Value from Spectronaut column 'EG.ProteinPTMLocations'

    Returns: List of target-mod positions only, for each protein ID sequence maps to.

    """
    sty_locations = []

    isoforms = PTM_ProteinLocations.split(';')
    for iso in isoforms:
        sites = iso.replace('(','').replace(')','').split(',')

        sty_only = []
        for s in sites:                         #Only keep position information from STY amino acids
            if s.startswith(('S','T','Y')):
                sty_only.append(s)

        sty_locations.append(sty_only)

    return(sty_locations)


def generate_collapsed_report(report):
    """
    Takes normal report from Spectronaut output and generates collapsed PTM report where each entry represents a unique phosphosite. A phosphosite
    is defined as a uniquely phosphorulated peptide sequence. Differentially heavy-modified peptides are listed as unique phosphosite entries.

    Args:
        report: Spectronaut Normal Report. Must contain the following columns for script to function:
            'R.FileName'
            'FG.IntMID'
            'EG.IntPIMID'
            'PG.ProteinAccessions'
            'EG.ProteinPTMLocations'
            'FG.Quantity'
            'EG.PTMAssayProbability'
            'EG.Cscore'

    Returns: Collapsed PTM site report. Additional columns added:
        'Only_Target_Mods': Peptide sequence with only target PTM annotated
        'HeavyMod_ProteinPositions': Protein location of heavy-labeled amino acids

    """

    # Only keep sequences that contain the target modification ([+80])
    target = report.loc[report['FG.IntMID'].str.contains('\[\+80\]')]

    target_mods = []
    heavy_ProteinLocations = []
    collapse_keys = []

    for index, row in target.iterrows():
        all_mods = row['FG.IntMID']                 #This value annotates all mods, including heavy mods
        non_heavy_mods = row['EG.IntPIMID']


        PEP_pos = row['PEP.PeptidePosition']

        target_mods_seq = remove_experimental_mods(all_mods)  #Remove M-oxidation, carbamidomethylation, and N-term acetylation, INCLUDE HEAVIES
        only_targets_seq = remove_experimental_mods(non_heavy_mods)

        target_mods.append(only_targets_seq)   #Add altered sequence to new column

        position_info = find_position(non_heavy_mods, target_mods_seq, PEP_pos)
        heavy_ProteinLocations.append(position_info)         #[L[+10]564]


        protein_ids = row['PG.ProteinAccessions']
        PTM_ProteinLocations = row['EG.ProteinPTMLocations']                                #Existing column also contains PTM locations M-Ox, Cam, etc
        STY_ProteinLocations = str(sty_ProteinLocations(PTM_ProteinLocations))              #Create column that only contains PTM locations of STY's, this way even sequences with other modifications will still collapse if they have the same phosphorylation

        k = (protein_ids,STY_ProteinLocations,str(position_info))         #Collapse key contains: Protein IDs, PTM protein locations, heavy mod protein locations
        collapse_keys.append(k)



    target['Only_Target_Mods'] = target_mods                               #Creates new column with sequences only consisting of target PTM annotations
    target['HeavyMod_ProteinPositions'] = heavy_ProteinLocations           #Creates new column with heavy-mod protein locations
    target['Collapse'] = collapse_keys                                     #Creates new column with collapse key values


    #In this section, each group under a collapse key is further split into members associated with each individual sample, then quant is separately summed per run within the gorup and added to a dictionary
    #Each collapse key is a key in a dictionary where the values

    groups = target.groupby('Collapse')

    meta = {}
    for name, group in groups:
        meta[name] = None

        file_quant = {}

        for s in samples:           #Sum quant by file
            file_quant[s] = None

            sub = group[group['R.FileName'] == s]       #Individual run/file/sample
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

    #Select columns to keep in the ouput
    keep = ['R.Condition','R.FileName', 'R.Replicate', 'PG.Genes','PG.Organisms', 'PG.ProteinAccessions', 'PG.IBAQ', 'PG.Quantity', 'PEP.PeptidePosition', 'PEP.Quantity', 'EG.IntPIMID', 'EG.ProteinPTMLocations','EG.PTMPositions [Phospho (STY)]',
            'EG.PTMAssayProbability','EG.PTMLocalizationProbabilities','EG.Cscore','EG.IntCorrScore','FG.Charge','FG.IntMID','FG.CalibratedMassAccuracy (PPM)', 'Only_Target_Mods', 'HeavyMod_ProteinPositions']

    #Append file names as column titles, these columns contain file-specific quant values
    for s in samples:
        keep.append(s)

    #Generate new collapsed report with selected columns only
    collapsed_report = collapsed_MaxCScore[keep]
    collapsed_report.to_csv(wd + 'VM_CollapsedReport.tsv', sep='\t', index=False)
    return(collapsed_report)

def summarize_any(full, collapsed_report, sample):
    """
    Generates summary stats for any data frame (could be subset of original report)

    Args:
        full: Data frame that has not been collapsed
        collapsed_report: Collapsed PTM site report
        sample: Column name that will be used for the summary document

    Returns: Dictionary with values for
        % precursors containing target modification (STY)
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique STY sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed phosphosites, including heavy mods

    """
    all_sequences = full['FG.IntMID']  # All precursors
    num_modified_sequences = len(all_sequences.unique())                           # All precursors, including experimentally-introduced mods and heavy mods

    phos_sequences = [x for x in all_sequences if '[+80]' in x]  # Only phospho-containing precursors
    enrichment = (len(set(phos_sequences)) / num_modified_sequences) * 100         #Phospho-containing precursors / all precursors

    # Number of unique STY modified sequences (Heavy mods allowed, no experimentally-introduced mods)
    phospho_no_experimental = list(map(remove_experimental_mods, phos_sequences))  # Remove experimentally-introduced modifications
    unique_phosphopeptides = len(set(phospho_no_experimental))                     #This will represent the number of phosphopeptides

    #For individual samples, find number of phosphosites
    if sample != 'Combined':
        select = collapsed_report[sample].values.tolist()
        res = [i for i in select if i is not None]
        phosphosites = len(res)

    #Length of data frame is the number of phosphosites identified across all samples searched together
    if sample == 'Combined':
        phosphosites = len(collapsed_report)

    return({'Run': sample, '%_Precursors_with_sty': enrichment, 'num_modified_sequences': num_modified_sequences ,'num_phosphopeptides': unique_phosphopeptides, 'num_phosphosites': phosphosites})

def summarize_non_missing(setlist, collapsed_report, samples, sample):
    """
    Reports summary stats sequences present across all samples searched together.

    Args:
        setlist: list of sets, one set of sequences per sample
        collapsed_report: Collapsed PTM site report
        samples: List of all names
        sample: Column name that will be used for the summary document

    Returns: Dictionary with values for
        % precursors containing target modification (STY)
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique STY sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed phosphosites, including heavy mods

    """
    intersection = set.intersection(*setlist)                               #Sequences that are found across all samples
    num_modified_sequences = len(intersection)

    phos_sequences = [x for x in intersection if '[+80]' in x]  # Only phospho-containing precursors
    enrichment = (len(set(phos_sequences)) / num_modified_sequences) * 100  #Phospho-containing precursors / all precursors

    phospho_no_experimental = list(map(remove_experimental_mods, phos_sequences))  #Remove experimentally-introduced modifications, keep heavy mods
    unique_phosphopeptides = len(set(phospho_no_experimental))  #This will represent the number of phosphopeptides

    no_missing = collapsed_report.dropna(subset = samples)
    no_missing.to_csv(wd + 'Intersection.tsv', sep = '\t', index = False)
    phosphosites = len(no_missing)

    return({'Run': sample, '%_Precursors_with_sty': enrichment, 'num_modified_sequences': num_modified_sequences,'num_phosphopeptides': unique_phosphopeptides, 'num_phosphosites': phosphosites})




def summarize(report, collapsed_report, samples):
    """
    Reports summary stats across all runs searched together.

    Args:
        report: Spectronaut Normal Report
        collapsed_report: Collapsed PTM site report
        samples: List of all file names

    Returns: Datafrane with following values for each sample, combined data set, and intersection across samples:
        % precursors containing target modification (STY)
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique STY sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed phosphosites, including heavy mods

    """
    print('Summarizing data...')
    data = pd.DataFrame(columns= ['Run','%_Precursors_with_sty', 'num_modified_sequences', 'num_phosphopeptides', 'num_phosphosites'])       #Initialize datarframe with columns

    non_missing = []                #Collect precursors found in each sample individually

    for s in samples:
        one = report[report['R.FileName'] == s]
        ret = summarize_any(one,collapsed_report, s)
        data = data.append(ret, ignore_index= True)

        set_IDs = set(one['FG.IntMID'])
        non_missing.append(set_IDs)

    #Add row for summarized combined data
    total = summarize_any(report, collapsed_report, 'Combined')
    data = data.append(total, ignore_index= True)

    #Add row for intersection of data across all samples searched together
    intersection = summarize_non_missing(non_missing, collapsed_report, samples, 'Intersection')
    data = data.append(intersection, ignore_index= True)

    data.to_csv(wd + 'VMSummary.tsv', sep='\t', index=False)
    return(data)


if __name__ == "__main__":
    wd = 'Y:\\LabMembers\\Tan\\CompGroup\\PTMSiteCollapse\\Pro\\'
    file_name = '20221209_145002_20221004_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_SmallLib0.75Loc_LocalizationScores_Report.tsv'

    info = prepare(wd, file_name)
    report = info[0]
    samples = info[1]

    collapsed_report = generate_collapsed_report(report)
    summarize(report, collapsed_report, samples)





