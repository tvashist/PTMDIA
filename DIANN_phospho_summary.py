import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os

#Select which instrument you're working on data from
# platform = 'Exploris'
# platform = 'Exploris_FAIMS'
platform = 'TimsTOF_Pro'
# platform = 'TimsTOF_SCP'


#Select Library
library = 'Library_3SS_Spiked'
# library = 'Library_12fxn_NotSpiked'
# library = 'Library_Combined'


report_directory_path = 'Z:/Helium_Tan/FINAL_PTMDIA/' + platform + '/DIANN/' +library +'/SearchOutputs/dia_nn/out/'     #Where is your Spectronaut output report?
diann = pd.read_csv(report_directory_path + 'report.tsv', delimiter= '\t', low_memory= False)


print("read")

#Get a summary of the number of phosphopeptide and phosphosite counts in the runs without spiked peptides
descriptor = platform + '_PhosphoBG_NoSpike'            #How you want your replicates to be annotated


nospike = diann.loc[diann['Run'].str.contains('NoSpike')]        #Only look at the entries that are from the no spike runs

files = list(set(nospike['Run']))                                                  #Unique runs in file, so all three reps of no spike runs
print(files)

phospho_summary = {}
phospho_summary['Run'] = []
phospho_summary['Phosphopeptides'] = []

pY = 0
for rep in range(0, len(files)):
    unique_phospho_status = []

    run = files[rep]

    phosphos = nospike.loc[nospike['Run'] == run]                            #Subset to just one rep
    phosphos = phosphos.loc[phosphos['Modified.Sequence'].str.contains('UniMod:21')]          #All entries that contain phosphorylated amino acids

    for x in phosphos['Modified.Sequence']:
        new = x.replace('(UniMod:1)','').replace('(UniMod:4)','').replace('(UniMod:35)','').replace('(UniMod:267)','').replace('(UniMod:259)','')                          #Remove annotations for oxidized methionine and carbamidomethylation
        unique_phospho_status.append(new)
        if 'Y(UniMod:21)' in new:
            pY += 1




    unique_phosphopeptides = len(set(unique_phospho_status))                        #Sequences with unique phosphorylation status
    print(pY)
    print(unique_phosphopeptides)




    phospho_summary['Run'].append(str(rep+1))
    phospho_summary['Phosphopeptides'].append(unique_phosphopeptides)
    # phospho_summary['Phosphosites'].append(phosphosites)


report = pd.DataFrame.from_dict(phospho_summary)

path = report_directory_path + platform + '_PhosphoStats/'
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)

report.to_csv(path + platform + '_PhosphoStatsSummary.tsv', sep='\t')







