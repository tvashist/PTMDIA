import pandas as pd
import re
import math
from statistics import mean,stdev,mode

#Get a summary of the number of phosphopeptide and phosphosite counts in the runs without spiked peptides
spectronaut = pd.read_csv('20220714_100348_PTMDIAProject_Pro_PhosphoBG_Report.tsv', delimiter='\t')
# head = spectronaut[:50]
# head.to_csv('ShortPhoshpoBGSCP.tsv', sep = '\t')

nospike = spectronaut.loc[spectronaut['R.FileName'].str.contains('NoSpike')]        #Only look at the entries that are from the no spike runs
reps = set(nospike['R.FileName'])

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

    return(len(set(sites)))

for run in reps:
    phosphos = nospike.loc[nospike['R.FileName'] == run]
    phosphos = phosphos.loc[phosphos['EG.IntPIMID'].str.contains('[+80]')]
    print(run[-4:])
    unique_phosphopeptides = len(set(phosphos['EG.IntPIMID']))
    phosphosites = count_sites(phosphos)
    print("Number of unique phosphopeptides: " + str(unique_phosphopeptides))
    print("Number of phosphosites: " + str(phosphosites))




