import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

#Uncomment which instrument you're working on data from#
# platform = 'Exploris'
# platform = 'Exploris_FAIMS'
# platform = 'TimsTOF_Pro'
platform = 'TimsTOF_SCP'

#Library#
library = 'Library_3SS_Spiked'
# library = 'Library_12fxn_NotSpiked'
# library = 'Library_Combined'
# library = 'directDIA'

# report_directory_path = 'Z:/Helium_Tan/FINAL_PTMDIA/' + platform + '/Spectronaut/' + library + '/SearchOutputs/'
report_directory_path = 'Z:/Alexi/TimsTOF_Jurkat_iRTs_phospho/'
report = '20221117_155033_Phospho_Jurkat_with_Pro12FxLibv1_Report_Version2_xls.tsv'
spectronaut = pd.read_csv(report_directory_path + report, delimiter='\t', low_memory=False)


print("read")
#Get a summary of the number of phosphopeptide and phosphosite counts in the runs without spiked peptides
# descriptor = platform + '_PhosphoBG_NoSpike'            #How you want your replicates to be annotated
descriptor = "Pro_12fxnLibrary"

# nospike = spectronaut.loc[spectronaut['R.FileName'].str.contains('NoSpike')]        #Only look at the entries that are from the no spike runs
nospike = spectronaut
reps = list(set(nospike['R.FileName']))                                                  #Unique runs in file, so all three reps of no spike runs
# nospike.to_csv(report_directory_path + "NoSpikesOnly.tsv", sep = '\t')
print(reps)
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
# phospho_summary['Phosphosites'] = []



for rep in range(0, len(reps)):
    unique_phospho_status = []

    run = reps[rep]
    # run_id = descriptor +"_"+ "0" + str(rep+1)                                        #Run descriptor that will go in data frame
    run_id = run
    phosphos = nospike.loc[nospike['R.FileName'] == run]                            #Subset to just one rep
    phosphos = phosphos.loc[phosphos['EG.IntPIMID'].str.contains('80')]          #All entries that contain phosphorylated amino acids
    for x in phosphos['EG.IntPIMID']:
        new = x.replace('[+16]','').replace('[+57]','').replace('[+42]','').replace('[+8]','').replace('[+10]','')                             #Remove annotations for oxidized methionine and carbamidomethylation
        unique_phospho_status.append(new)
        # print(new)
    # print("I just changed all the mod annotations")


    unique_phosphopeptides = len(set(unique_phospho_status))                        #Sequences with unique phosphorylation status

    # phosphosites = count_sites(phosphos)                                            #Unique phosphosites based on location on protein sequence

    phospho_summary['Run'].append(run_id)
    phospho_summary['Phosphopeptides'].append(unique_phosphopeptides)
    # phospho_summary['Phosphosites'].append(phosphosites)


report = pd.DataFrame.from_dict(phospho_summary)

path = report_directory_path + platform + '_PhosphoStats/'
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)

report.to_csv(path + platform + descriptor + '_PhosphoStatsSummary.tsv', sep='\t')







