import pandas as pd
import re
import math
from statistics import mean,stdev,mode

#Get a summary of the number of phosphopeptide and phosphosite counts in the runs without spiked peptides
descriptor = 'Pro_PhosphoBG_NoSpike'            #How you want your replicates to be annotated
spectronaut = pd.read_csv('C:/Users/tvashist/PycharmProjects/PTMDIA_Project/20220714_100348_PTMDIAProject_Pro_PhosphoBG_Report.tsv', delimiter='\t')     #Entire data frame

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


print("here before loop!")
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
report.to_csv('C:/Users/tvashist/PycharmProjects/PTMDIA_Project/PhosphoBG_Curve/Pro_NoSpike_Summary/Summary.tsv', sep = '\t')





