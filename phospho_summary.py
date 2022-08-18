import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

#Select which instrument you're working on data from
# # platform = 'Exploris'
# platform = 'Pro'
platform = 'SCP'


report_directory_path = 'Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/' + platform + "/"     #Where is your Spectronaut output report?

if platform == 'SCP':
    spectronaut = pd.read_csv(report_directory_path + '20220818_155021_PTMDIAProject_DIACurveAnalysis_WithSpecLib_Report.tsv', delimiter='\t')

if platform == 'Pro':
    spectronaut = pd.read_csv(report_directory_path + '20220803_134456_PTMDIAProject_DIACurveAnalysis_WithSpecLib_Report_addedFGLabel.tsv', delimiter= '\t', low_memory= False)

if platform == 'Exploris':
    spectronaut = pd.read_csv('', delimiter = '\t')

print("read")
#Get a summary of the number of phosphopeptide and phosphosite counts in the runs without spiked peptides
descriptor = platform + '_PhosphoBG_NoSpike'            #How you want your replicates to be annotated

nospike = spectronaut.loc[spectronaut['R.FileName'].str.contains('NoSpike')]        #Only look at the entries that are from the no spike runs
reps = list(set(nospike['R.FileName']))                                                  #Unique runs in file, so all three reps of no spike runs


def count_sites(df):

    sites = []
    for index, row in df.iterrows():
        protein = row['PG.ProteinAccessions'].split(';')[0]                                     #Only select first protein
        if type(row['EG.ProteinPTMLocations']) != float:
            site = row['EG.ProteinPTMLocations'].split(';')[0].replace('(','').replace(')','')      #Site of the phosphorylation in the PROTEIN, not peptide
        site_new = site.split(',')                                                              #If protein contains multiple sites, this will make a list of them

        for x in range(len(site_new)):
            site_info = (protein, site_new[x])              #Append info of protein and the site such that each value is a unique phosphosite

            if site_info not in sites:
                sites.append(site_info)

    return(len(set(sites)))                                 #Returns number of unique phosphosites based on location on protein sequence

phospho_summary = {}
phospho_summary['Run'] = []
phospho_summary['Phosphopeptides'] = []
phospho_summary['Phosphosites'] = []



for rep in range(0, len(reps)):
    unique_phospho_status = []

    run = reps[rep]
    run_id = descriptor +"_"+ "0" + str(rep)                                        #Run descriptor that will go in data frame

    phosphos = nospike.loc[nospike['R.FileName'] == run]                            #Subset to just one rep
    phosphos = phosphos.loc[phosphos['EG.IntPIMID'].str.contains('[+80]')]          #All entries that contain phosphorylated amino acids
    for x in phosphos['EG.IntPIMID']:
        new = x.replace('[+16]','').replace('[+57]','')                             #Remove annotations for oxidized methionine and carbamidomethylation
        unique_phospho_status.append(new)
    print("I just changed all the mod annotations")

    unique_phosphopeptides = len(set(unique_phospho_status))                        #Sequences with unique phosphorylation status
    phosphosites = count_sites(phosphos)                                            #Unique phosphosites based on location on protein sequence

    phospho_summary['Run'].append(run_id)
    phospho_summary['Phosphopeptides'].append(unique_phosphopeptides)
    phospho_summary['Phosphosites'].append(phosphosites)


report = pd.DataFrame.from_dict(phospho_summary)

path = report_directory_path + platform + '_PhosphoStats/'
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)

report.to_csv(path + platform + '_PhosphoStatsSummary.tsv', sep='\t')







