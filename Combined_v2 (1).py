#!/usr/bin/env python
# coding: utf-8

# In[4]:


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[78]:


from alphamap.importing import import_data
from alphamap.preprocessing import format_input_data
from pyteomics import fasta

import re
import math
import os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

from statistics import mean,stdev,mode
import argparse
import sys

import yaml
from yaml.loader import SafeLoader


# In[112]:


wd = 'Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Combined_DDA_Predicted_Library\dia_nn\out' + '\\'
with open(wd + 'input.yaml', 'r') as file:
    params = yaml.safe_load(file)


# In[113]:


#Unpack yaml file input parameters
file_name = params['main']['report']
fasta_path = params['main']['fasta_path']   #File path that points directly to the file
engine = params['main']['engine']

target_ptm = params[engine]['target_ptm']
experimental = params[engine]['experimental']
heavies = params[engine]['heavies'] 
file_col = params[engine]['file']
all_sequences = params[engine]['all_sequences']
prot_column = params[engine]['prot_column']

ptm_phrase = params['phrases'][target_ptm]
threshold = params[engine]['ptm_confidence_threshold']

print(threshold, engine)


# In[9]:


#DIANN Pre-processing#
def preprocess(report):
    """
    Localizes report for site confidence >= threshold, exports localized report. Performs alphamap processing of report to get peptide position and PTM information. Adjusts PTM annotations.

    Args: Long-form DIA-NN report

    Returns: (Localized dataframe report, alphamap output) for further processing of DIA-NN data
    """
    #Localize report, export localized report
    if threshold != 0:
        localized_report = report[report['PTM.Site.Confidence'] >= threshold]
        localized_report.to_csv(wd + 'Localized_Report.tsv', sep = '\t', index = False)
        report_path = wd + 'Localized_Report.tsv'
        
    else: 
        localized_report = report
        report_path = wd + file_name
        

    print('Report with confidently localized entries exported.')

    #Prepare report and fasta for Alphamap processing
    diann_report = import_data(report_path)
    fasta_file = fasta.IndexedUniProt(fasta_path)
    print('Prepared for alphamap processing')

    #Remove contaminants and heavy annotations (does not remove full sequences containing heavies)
    diann_report = diann_report[diann_report["all_protein_ids"].str.contains("contaminant") == False].reset_index(drop=True)
    diann_report = diann_report.replace('\(UniMod:259\)','', regex = True)
    diann_report = diann_report.replace('\(UniMod:267\)','', regex = True)

    #Run report throught alphamap to get formatted data with peptide locations
    formatted_data = format_input_data(df=diann_report, fasta=fasta_file, modification_exp = r'\[.*?\]')  # Run through alphamap to get peptide locations on protein
    formatted_data.rename(columns={'modified_sequence': 'Modified.Sequence', 'all_protein_ids': "Protein.Ids"}, inplace=True)   #Rename 'all_protein_ids' column to 'Protein IDs to match the matrix it will be merged with

    #Update peptide start and end positions to be 1-indexed instead of 0-indexed
    formatted_data['start'] = formatted_data['start'] + 1
    formatted_data['end'] = formatted_data['end'] + 1

    #Change annotations in the formatted data to match matrix it will be merged with
    formatted_data = formatted_data.replace('\[Phospho \(STY\)\]','(UniMod:21)', regex=True)
    formatted_data = formatted_data.replace('\[Carbamidomethyl \(C\)\]','(UniMod:4)', regex=True)
    formatted_data = formatted_data.replace('\[Acetyl \(N-term\)\]','(UniMod:1)', regex=True)
    formatted_data = formatted_data.replace('\[Oxidation \(M\)\]','(UniMod:35)', regex=True)

    formatted_data = formatted_data.sort_values('unique_protein_id')        #Sort unique protein IDs so when STY locations get appended, they match order of alphabetized protein IDs
    
    return(localized_report, formatted_data)


# In[10]:


def filter_matrix(localized_report, matrix): #DIA-NN
    """
    Filters DIA-NN matrix file-specific quant for sequences that are only in the localized report.
    
    Args: Localized report, pre-collapsed matrix
    
    Returns: Matrix where intensities for each sample are only reported on file-specific localized sequences.
    
    """
    
    localized = {}
    for s in samples:           
        one = localized_report[localized_report['File.Name'] == s]   #Only entries from that run in report
        if s not in localized:                   #Key is sample name, list of unique modified sequences that are localized in that run
            localized[s] = None
            localized[s] = one['Modified.Sequence'].unique()
            
        else:
            localized[s] = one['Modified.Sequence'].unique()
            
    print('Found localized')
            
    for index, row in matrix.iterrows():          #Only iterate through matrix once
        seq = row['Modified.Sequence']
        
        for s in samples:         
            if seq not in localized[s]:
                matrix.at[index,s] = None         #If the sequence in the matrix is not localized, replace intensity with None
                
    return(matrix)


# In[114]:


if engine == 'SN':
    report = pd.read_csv(wd + file_name, delimiter ='\t', low_memory= False)        #Read in Spectronaut "Normal Report"
    if threshold != 0:
        report = report.loc[report['EG.PTMAssayProbability'] >= threshold]


    samples = (report[file_col]).unique()                                       #List of unique samples/ LC-MS runs that IDs may come from

    
elif engine == 'DIANN':
    report = pd.read_csv(wd+ file_name, sep = '\t')
    matrix = pd.read_csv(wd + 'report.pr_matrix.tsv', sep = '\t', low_memory = False)
    samples = (report[file_col]).unique()
    
    print('Running DIA-NN report through alphamap library')
    processed = preprocess(report)
    localized_report = processed[0]
    formatted_data = processed[1]
    
    filtered_matrix = filter_matrix(localized_report, matrix)
    
else:
    raise Exception("Sorry, please check the spelling of the software.")
    
print('Ready for processing.')


# In[12]:


def remove_experimental_mods(sequence):   #SN and DIA-NN
    """
    Removes annotations relating to N-term acetylation, oxidation of methionine, and carbamidomethylation. Also removes '_' flanking sequences to ease indexing in Spectronaut.

    Args:
        sequence: Sequence containing experimentally-introduced mod annotations listed above.

    Returns: Peptide sequence without experimentally-introduced mod annotations listed above.
    """
    for item in experimental:
        sequence = sequence.replace(item, '') #Remove carbamidomethylation, N-term acetylation, M-oxidation
        
    return(sequence)


# remove_experimental_mods('_EGPTDHLESAC[+57]PL[+8]NLPLQNNHTAADM[+16]YLS[+80]PVR[+10]S[+80]PK_')


# In[13]:


def generate_collapsed_report_SN(report):
    """
    Takes normal report from Spectronaut output and generates collapsed PTM report where each entry represents a unique phosphosite. A phosphosite
    is defined as a uniquely phosphorylated peptide sequence. Differentially heavy-modified peptides are listed as unique phosphosite entries.

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

    Returns: Collapsed PTM site report.
    """

    target = report.loc[report['EG.ModifiedSequence'].str.contains(ptm_phrase)]    #Only entries with target PTM are filtered

    
    collapse_keys = []

    targets_unannotated = []
    for index, row in target.iterrows():
        
        sequence = remove_experimental_mods(row['FG.IntMID'])  
        pep_prot_positions = row['PEP.PeptidePosition']
        protein_ids = row['PG.ProteinAccessions']
        
        #Target and heavy PTM positions
        ptm_prot_positions = find_prot_positions(sequence, pep_prot_positions)
        
        #Integer-unannotated target positions used by find_flanking function
        unannotated = [[re.sub("\[.*?\]", "", s) for s in pos if target_ptm in s] for pos in ptm_prot_positions]  #Keeps only sequences with target mod (ex: '+80'), and removes integer annotations in brackets
        targets_unannotated.append(unannotated)
        
        #Collapse keys
        k = (protein_ids,str(ptm_prot_positions))    #Collapse key contains: Protein IDs, target PTM protein locations (this way even sequences with other modifications will still collapse if they have the same target PTM positions), heavy mod protein locations
        collapse_keys.append(k)

        
    #Add columns to target dataframe
    target['Collapse'] = collapse_keys                                     #Creates new column with collapse key values
    target[ptm_phrase +'ProteinLocations'] = targets_unannotated           #Creates new column with integer-unannotated target PTM protein positions
    
    first_column = target.pop('Collapse')
    target.insert(0, 'Collapse', first_column)
    
#     target.to_csv(wd + 'Keyed_PreCollapse.tsv', sep = '\t', index = False)
    
    #In this section, each group under a collapse key is further split into members associated with each individual sample, then quant is separately summed per run within the gorup and added to a dictionary
    #Each collapse key is a key in a dictionary where the values

    groups = target.groupby('Collapse')   #Group entries by collapse key
    
    #Structure: {Collapse_Key: {Sample1: Summed_quant_01, Sample2: Summed_quant_02}}


    meta = {}
    for n, group in groups:
        
        meta[n] = None              #n is collapse key

        file_quant = {}

        for s in samples:           #Sum quant by file
            file_quant[s] = None

            sub = group[group['R.FileName'] == s]       #All entries from that individual run/file/sample
            if len(sub) > 0:
                quant = sub['FG.Quantity'].sum()              #Sum quant
            else:
                quant = None                                  #Otherwise append None

            file_quant[s] = quant

        meta[n] = file_quant                                  #Adding quant dict for that collapse key

    #In the following section, rows are being collapsed
    
    #Only keep row from a group that contains the highest PTM Assay Probability Score (most confident localization); There will be repeats for rows with same PTM Assay Probability Score
    group_MaxPTMAssayScore = target.groupby(['Collapse'])['EG.PTMAssayProbability'].transform(max) == target['EG.PTMAssayProbability']
    collapsed_MaxPTMAssayScore = target[group_MaxPTMAssayScore]

    #Make collapse key first column, will remove this helper column for final output
    first_column = collapsed_MaxPTMAssayScore.pop('Collapse')
    collapsed_MaxPTMAssayScore.insert(0, 'Collapse', first_column)

    #Only keep row from group with highest EG.Cscore; There will be no more repeats
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
    keep = ['Collapse','R.Condition','R.FileName', 'R.Replicate', 'PG.Genes','PG.Organisms', 'PG.ProteinAccessions', 'PG.IBAQ', 'PG.Quantity', 'PEP.PeptidePosition', ptm_phrase +'ProteinLocations', 'PEP.Quantity', 'EG.IntPIMID', 'EG.ProteinPTMLocations',
            'EG.PTMAssayProbability','EG.PTMLocalizationProbabilities','EG.Cscore','EG.IntCorrScore','FG.Charge','FG.IntMID','FG.CalibratedMassAccuracy (PPM)']

    #Append file ptm_phrases as column titles, these columns contain file-specific quant values
    keep.extend(samples)

    #Generate new collapsed report with selected columns only
    collapsed_report = collapsed_MaxCScore[keep]
    collapsed_report.to_csv(wd + 'VM_CollapsedReport.tsv', sep='\t', index=False)

    return(collapsed_report)


# In[14]:


def find_prot_positions(sequence,prot_positions):
    
    all = list(heavies)
    all.append(target_ptm)

    if engine == 'SN':
        all_regex = '[A-Z]('+ '|'.join(['\\[\\' + x + '\\]' for x in all]) + ')'        #Generates regex expression that searches for all target ptm and heavy annotations and the amino acid
        prot_positions = re.split(r';|,', prot_positions)
    
    if engine == 'DIANN':
        all_regex = '[A-Z]('+ '|'.join(['\\(' + x + '\\)' for x in all]) + ')'

    finder = re.finditer(all_regex, sequence)                  #Find matches within sequence
    naked = re.sub(all_regex, 'X', sequence)                   #Replace all matches with a single character to find their position in naked sequence

    annotations = [match.group() for match in finder]                               #All annotations with amino acid residue
    pep_positions = [match.start() for match in re.finditer('X', naked)]            #Peptide positions of matches (0-indexed)
        
    full = []
    
    for p in prot_positions:                                             #PEP_positions relating to different protein IDs are separated by ';'
        
        protein_pos = [str(int(p) + x) for x in pep_positions]                      #Add the peptide position of the match to the protein position of peptide --> protein position of the match
        ann = [''.join(item) for item in zip(*[annotations, protein_pos])]          #Combine the annotation (regex match) with the updated protein position
        full.append(ann)
        
    return(full)
        
# find_prot_positions('AAVSHWQQQS(UniMod:21)YLDSGIHSGATTTAPSLSGK(UniMod:259)', [40])
# find_prot_positions('TAGTSFMMT[+80]PYVVT[+80]R[+10]', '175;68;175;175')


# In[15]:


def condense_formatted_data(formatted_data):  #DIA-NN   
    """
    Adds PTM residue protein positions and condenses alphamap output such that PTM protein locations are given for all protein IDs a sequence maps to in a single row.
    
    Args: Alphamap output where each entry corresponds to a unique protein ID a sequence maps to.
    
    Returns: Formatted output collapsed with PTM protein locations for all protein IDs a sequence maps to.
    
    """
    
    #The Alphamap report creates individual entry for every protein ID a sequence is mapped to, here we group by 'Protein.ids' (all protein IDs a sequence maps to) and 'Modified Sequence'
    condense = ['start','end']
    list_cols = {col: lambda x: list(x) for col in condense}               #Condenses entries from individual rows to a list of PTM protein locations
    
    keep = ['Protein.Ids','Modified.Sequence','naked_sequence','PTMtypes','PTMsites']
    other_cols = {col: 'first' for col in keep}                                                           #Keeps first of all other entries because they are identical across all rows of the same group

    col_map = {**other_cols, **list_cols}                                                                 #Determines how grouping should be done


    formatted_data = formatted_data.groupby(['Protein.Ids','Modified.Sequence']).agg(col_map).reset_index(drop=True)
    
    return(formatted_data)


# In[16]:


def alphabetize(proteins):  #DIA-NN
    """
    Alphabetizes order of protein IDs in DIA-NN matrix 'Protein.Ids' column. This is needed for the merging of the alphamap output to the matrix.

    Args: string of protein IDs separated by ';'

    Returns: string of alphabetized protein IDs separated by ';'
    """
    ls = sorted(proteins.split(';'))    #Split and sort
    str = ';'.join(ls)                  #Rejoin with semicolon into str

    return(str)


# In[17]:


def merge(formatted_data, filtered_matrix):  #DIA-NN
    """
    Merges alphamap output with protein positions to the filtered intensity matrix, which does not have protein positions. 
    
    Args: Formatted alphamap output, DIA-NN matrix output
    
    Returns: Merged ouput of intensity matrix with protein position information appended.

    """
    print('Merging alphamap output and DIA_NN Precursor matrix')

    filtered_matrix['Precursors_Heavies'] = filtered_matrix['Modified.Sequence']       #Maintains record or heavies included in sequence
    
    for h in heavies:
        h = '(' + h + ')'
        filtered_matrix['Modified.Sequence'] = filtered_matrix.apply(lambda x: x['Modified.Sequence'].replace(h,''),axis=1) #Removes heavy annotations from this sequence so that it will match the sequences in alphamap output (which cannot contain heavy annotations)
        
    filtered_matrix['Protein.Ids'] = filtered_matrix.apply(lambda x: alphabetize(x['Protein.Ids']),axis=1)  # Alphabetizes protein IDs in matrix so
    

    #Make changes to alphamap output
    formatted_data = condense_formatted_data(formatted_data)                           #Adds PTM protein locations to formatted data

    new_df = filtered_matrix.merge(formatted_data, "inner", on=["Modified.Sequence", "Protein.Ids"])     #Merges two data frames based on common protein IDs and modified sequences that do not contain heavies
    new_df.to_csv(wd + 'Merged.tsv', sep = '\t')

    return(new_df)


# In[18]:


def generate_collapsed_report_DIANN(formatted_data, filtered_matrix):
    """
    Generates collapsed PTM report where each entry represents a unique phosphosite. A phosphosite
    is defined as a uniquely phosphorylated peptide sequence. Differentially heavy-modified peptides are listed as unique phosphosite entries.

    Args:
        formatted_data = alphamap output
        filtered_matrix = Filtered DIA_NN precursor matrix


    Returns: Collapsed PTM site report. Additional columns added:
        'Only_Target_Mods': Peptide sequence with only target PTM annotated
        'HeavyMod_ProteinPositions': Protein location of heavy-labeled amino acids

    """
    report = merge(formatted_data, filtered_matrix)                                      #Generate merged dataframe that will then be collapsed

    target = report.loc[report['Modified.Sequence'].str.contains(target_ptm)]            #Only phospho
    

    collapse_keys = []
    targets_unannotated = []

    for index, row in target.iterrows():
        protein_ids = row['Protein.Ids']
        
        sequence = remove_experimental_mods(row['Precursors_Heavies'])                    #This column contains sequences containing heavy mods
        pep_prot_positions = row['start']
        
        #Target and heavy PTM positions
        ptm_prot_positions = find_prot_positions(sequence, pep_prot_positions)
        
        #Integer-unannotated target positions used by find_flanking function
        unannotated = [[re.sub("\(.*?\)", "", s) for s in pos if target_ptm in s] for pos in ptm_prot_positions]  #Keeps only sequences with target mod (ex: '+80'), and removes integer annotations in brackets
        targets_unannotated.append(unannotated)

        #Generate collapse key: #Collapse key contains: Protein IDs, PTM protein locations, heavy mod protein locations
        k = (protein_ids,str(ptm_prot_positions))         
        collapse_keys.append(k)

    target['Collapse'] = collapse_keys                                     #Creates new column with collapse key values
    target[ptm_phrase +'ProteinLocations'] = targets_unannotated

    #Prepare columns for grouping
    keep = ['Collapse', 'Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description',
            'Proteotypic', 'Stripped.Sequence', 'Modified.Sequence', 'Precursors_Heavies','Precursor.Id', 'start', 'end',
            'PTMsites', 'PTMtypes', ptm_phrase+'ProteinLocations']

    other_cols = {col: 'first' for col in keep}
    quants = {col: lambda x: x.sum(skipna=False) for col in samples}       #Sum all the values that are not NA

    col_map = {**other_cols, **quants}

    collapsed = target.groupby(['Collapse'], as_index=False).agg(col_map)  #Group by collapse key, sum values by sample


    # #Append file names as column titles, these columns contain file-specific quant values
    keep.extend(samples)

    collapsed_report = collapsed[keep]
#     collapsed_report.to_csv(wd + 'VM_CollapsedReport.tsv', sep = '\t', index = False )
    
    return(report, collapsed_report)  #Return merged matrix and collapsed report


# In[19]:


def generate_collapsed_report():
    if engine == 'SN':
        collapsed_report = generate_collapsed_report_SN(report)
        
    if engine == 'DIANN':
        collapsed_report = generate_collapsed_report_DIANN(formatted_data, filtered_matrix)[1]
        
    return(collapsed_report)


# In[35]:


def summarize_non_missing(setlist, collapsed_report, samples, sample, threshold):
    """
    Reports summary stats of sequences present across all samples searched together.

    Args:
        setlist: list of sets, one set of sequences per sample
        collapsed_report: Collapsed PTM site report
        samples: List of all names
        sample: Column name that will be used for the summary document

    Returns: Dictionary with values for
        % precursors containing target modification 
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique PTM sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed PTM-sites, including heavy mods

    """
    intersection = set.intersection(*setlist)                                 #Sequences that are found across all samples
    num_modified_sequences = len(intersection)

    ptm_sequences = [x for x in intersection if target_ptm in x]              #Only PTM-containing precursors
    enrichment = (len(set(ptm_sequences)) / num_modified_sequences) * 100     #PTM-containing precursors / all precursors

    ptm_no_experimental = list(map(remove_experimental_mods, ptm_sequences))  #Remove experimentally-introduced modifications, keep heavy mods
    unique_ptm_peptides = len(set(ptm_no_experimental))                       #This will represent the number of phosphopeptides

    no_missing = collapsed_report.dropna(thresh = threshold, subset = samples)                    #Drop sites that have missing values across ANY of the samples 
    ptm_sites = len(no_missing)

    
    
    ret_df = pd.DataFrame({'Run': sample, '%_Precursors_with_' + ptm_phrase: enrichment,
             'num_modified_sequences': num_modified_sequences,
             'num_' + ptm_phrase + '_peptides': unique_ptm_peptides,
             'num_' + ptm_phrase + '_sites': ptm_sites}, index=[0])

    return (ret_df)


# In[34]:


def summarize_any(full, collapsed_report, sample):
    """
    Generates summary stats for any data frame (could be subset of original report)

    Args:
        full: Data frame that has not been collapsed
        collapsed_report: Collapsed PTM site report
        sample: Column name that will be used for the summary document

    Returns: Dictionary with values for
        % precursors containing target modification
        Number of modified peptides: Including different versions of experimentally-introduced modifications and heavy mods
        Number of phosphopeptides: Number of unique PTM sequences, NOT including experimentally-introduced modifications, but including heavy mods
        Number of phosphosites: Collapsed PTM-sites, including heavy mods

    """

    num_modified_sequences = len(full[all_sequences].unique())                          #All precursors, including experimentally-introduced mods and heavy mods
    
    #Enrichment efficiency
    ptm_sequences = [x for x in full[all_sequences] if target_ptm in x]                 #Only PTM-containing precursors
    enrichment = (len(set(ptm_sequences)) / num_modified_sequences) * 100         #PTM-containing precursors / all precursors
    

    #Number of unique PTM modified sequences
    ptm_no_experimental = list(map(remove_experimental_mods, ptm_sequences))      #Remove experimentally-introduced modifications, keep heavy mods
    unique_ptm_peptides = len(set(ptm_no_experimental))                           #This will represent the number of PTM-peptides

    #Number of phosphosites per sample
    if sample != 'Combined':
        select = collapsed_report[sample].values.tolist()
        
        if engine == 'SN':
            res = [i for i in select if i is not None]
            
        if engine == 'DIANN':
            res = [i for i in select if math.isnan(i) == False]
        ptm_sites = len(res)

    #Length of data frame is the number of phosphosites identified across all samples searched together
    if sample == 'Combined':
        ptm_sites = len(collapsed_report)
        
    #Create data frame to add to summary report
    ret_df = pd.DataFrame({'Run': sample, '%_Precursors_with_' + ptm_phrase: enrichment,
        'num_modified_sequences': num_modified_sequences,
        'num_' + ptm_phrase + '_peptides': unique_ptm_peptides,
        'num_' + ptm_phrase + '_sites': ptm_sites}, index=[0])
    
    
    return (ret_df)


# In[33]:


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
    data = pd.DataFrame(columns= ['Run', '%_Precursors_with_' + ptm_phrase, 'num_modified_sequences','num_' + ptm_phrase + '_peptides','num_' + ptm_phrase + '_sites'])       #Initialize datarframe with columns

    non_missing = []                                        #Collect precursors found in each sample individually, will be list of sets

    for s in samples:
        one = report[report[file_col] == s]
        ret = summarize_any(one, collapsed_report, s)       #Returns a row of the dataframe containing stats for that sample
        data = pd.concat([data, ret], ignore_index=True)

        set_IDs = set(one[all_sequences])
        non_missing.append(set_IDs)

    #Add row for summarized combined data
    total = summarize_any(report, collapsed_report, 'Combined')
    data = pd.concat([data, total], ignore_index=True)

    #Add row for intersection of data across all samples searched together
    intersection = summarize_non_missing(non_missing, collapsed_report, samples, 'Intersection',len(samples))
    fifty_perc_complete = summarize_non_missing(non_missing, collapsed_report, samples, 'Fifty Percent Complete',len(samples)/2)
    
    data = pd.concat([data, intersection], ignore_index=True)
    data = pd.concat([data, fifty_perc_complete], ignore_index = True)

    data.to_csv(wd + 'VMSummary_' + str(threshold) + '.tsv', sep='\t', index=False)
    print('Summarizing complete')
    
    return(data)


# In[29]:


def get_summary():
    
    if engine == 'SN':
        summarize(report, collapsed_report, samples)
        
    if engine == 'DIANN':
        summarize(localized_report, collapsed_report, samples)


# In[24]:


def find_flanking(collapsed_report, organism):  #Adapted for both software
    
    print('Finding Flanks!!!')
    
    meta = pd.read_csv('Z:\\Helium_Tan\\PhosphositePlus\\Phosphorylation_site_dataset', sep = '\t')
    meta = meta[meta['ORGANISM'] == organism]   #Filter for organism, typically human

    dict = {}

    for index, row in meta.iterrows():
        prot_id = row['ACC_ID']
        residue = row['MOD_RSD']
        flank_seq = row['SITE_+/-7_AA'].upper()

        res = residue.split('-')[0]
        k = (prot_id, res)

        if k not in dict:
            dict[k] = None
            dict[k] = flank_seq

        else:
            dict[k] = flank_seq



    print('Done building dictionary. Now appending flanked sequences to collapsed report...')


    flanks_column = []         #Global list
    
    #Iterate through entries and add flanking sequences for each protein ID and site
    for index, row in collapsed_report.iterrows():
        
        prot_id = row[prot_column].split(';')                             #Get all protein IDs a sequence is mapped to
        locations = list(row[ptm_phrase +'ProteinLocations'])             #PTM locations in protein

        row_flanks = {}                                             #Each row/entry will have a dictionary if there are flanking sequences in phosphosite table, with protein ID as key and flanking sequence(s) as value(s)
        
        #Find flanks for each protein ID listed in entry
        for i in range(0, len(prot_id)):                            
            
            prot_flanks = []         #If there are multiple sites on the same sequence, they will all be appended here
            
            sites = locations[i]     #Can contain multiple sites
            id = prot_id[i]

            contains = False
            for s in sites:                                          #Each PTM-site on sequence associated with that protein ID
                rep_k = (id,s)                                       #(Protein ID, PTM Location on that protein ID)

                if rep_k in dict.keys():                             #Phosphosite table has mapping seq
                    contains = True
                    prot_flanks.append(dict[rep_k])
                    
                else:
                    prot_flanks.append('')                           #Sometimes if a seq is doubly phosphorylated only one site will be mapped in phosphosite table


            if contains:
                if id not in row_flanks:
                    row_flanks[id] = None
                    row_flanks[id] = prot_flanks
                else:
                    row_flanks[id] = prot_flanks

        if len(row_flanks) > 0:
            flanks_column.append(row_flanks)
        else:
            flanks_column.append(None)


    collapsed_report['7AA_Flanking'] = flanks_column
    collapsed_report.to_csv(wd + 'VM_CollapsedReport_Flanks.tsv', sep = '\t', index = False)
    print('Completed adding flanking sequences.')
    
    return(collapsed_report)


# In[115]:


collapsed_report = generate_collapsed_report()
get_summary()

# if ptm_phrase == 'Phospho':
#     find_flanking(collapsed_report, "human")

