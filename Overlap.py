#!/usr/bin/env python
# coding: utf-8

# In[1]:


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[2]:


import pandas as pd
import matplotlib
from matplotlib_venn import venn2, venn3
from upsetplot import plot
import seaborn as sns
import numpy as np
import math


# In[3]:


#Spectronaut
SN_DDALib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\DDALibrary\\20230308_091202_R02_PhosphoDIA_Pro_DDALibrary_Report.tsv', sep = '\t', low_memory = False)
SN_directDIA = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\directDIA\\20230227_104026_PhosphoDIA_R02_Report_directDIA.tsv', sep = '\t', low_memory = False)
SN_HybridLib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\Hybrid\\20230313_103651_R02_PhosphoDIA_Pro_HybridLibrary_Report.tsv', sep = '\t', low_memory = False)


# In[19]:


DIANN_DDALib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\DDA_Library\Three_varmods\dia_nn\out\\report.tsv', sep = '\t')
DIANN_CombinedLib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Combined_DDA_Predicted_Library\dia_nn\out\\report.tsv', sep = '\t')
DIANN_PredictedLib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Predicted_Library\dia_nn\out\\report.tsv', sep = '\t')


# In[4]:


def convert(sequence):
    
    converted = sequence.replace('[+80]','(UniMod:21)').replace('_[+42]','(UniMod:1)').replace('[+57]','(UniMod:4)').replace('[+16]','(UniMod:35)').replace('_','').replace('[+8]','(UniMod:259)').replace('[+10]','(UniMod:267)')
    return(converted)


# In[5]:


SN = [SN_DDALib,SN_directDIA, SN_HybridLib]
for search in SN:
    search['Modified.Sequence'] = search['FG.IntMID'].apply(lambda x: convert(x))


# In[6]:


def filter_sn(df):
    phospho = df.loc[df['Modified.Sequence'].str.contains('\(UniMod:21\)')]
    filtered = phospho.loc[phospho['EG.PTMAssayProbability'] >= 0.99]
    return(filtered)


# In[7]:


SN_DDALib = filter_sn(SN_DDALib)
SN_directDIA = filter_sn(SN_directDIA)
SN_HybridLib = filter_sn(SN_HybridLib)


# In[20]:


def filter_diann(df):
    phospho = df.loc[df['Modified.Sequence'].str.contains('\(UniMod:21\)')]
    filtered = phospho.loc[phospho['PTM.Site.Confidence'] >= 0.51]
    print(len(df), len(phospho),len(filtered))                                                                
    return(filtered)


# In[21]:


DIANN_DDALib = filter_diann(DIANN_DDALib)
DIANN_CombinedLib = filter_diann(DIANN_CombinedLib)
DIANN_PredictedLib = filter_diann(DIANN_PredictedLib)


# In[22]:


# all = [SN_DDALib,SN_directDIA, SN_HybridLib]
all = [DIANN_DDALib, DIANN_PredictedLib, DIANN_CombinedLib]
# all = [DIANN_DDALib, DIANN_PredictedLib, DIANN_CombinedLib,SN_DDALib,SN_directDIA, SN_HybridLib]

sequences_all = [set(x['Modified.Sequence']) for x in all]


# In[23]:


# names = ['SN_DDALib','SN_directDIA', 'SN_HybridLib']
names = ['DIANN_DDALib','DIANN_PredictedLib', 'DIANN_HybridLib']
# names = ['DIANN_DDALib','DIANN_PredictedLib', 'DIANN_HybridLib','SN_DDALib','SN_directDIA', 'SN_HybridLib']


# In[24]:


venn3(sequences_all, names,set_colors=('#0000ff', '#ff4500','#006400'))


# In[167]:


def find_unique(all, sequences_all):
    
    updated = []
    
    for i in range (0, len(sequences_all)):        #For each dataframe
        seqs = sequences_all[i]
        
        new = list(sequences_all)  #Make copy of full list
        new.remove(seqs)

        other = set.union(*new)
        unique = seqs - other
        print(len(unique))

        #keep only sequences that are unique in the overarching dataframe
        
        df = all[i]
        find_unique = df.loc[df['Modified.Sequence'].isin(unique)]
        
        updated.append(find_unique)
    
    return(updated)
    
       

unique = find_unique(all, sequences_all)


# In[172]:


dict = {}
for i in range(0, len(all)):
    name = names[i]
    df = all[i]
    


    if 'SN' in name:
        df = df.dropna(subset=['EG.PTMAssayProbability'])

        
        charge = df['FG.Charge']
        im = df['FG.IonMobility']
        localization =df['EG.PTMAssayProbability']
        score = df['EG.Cscore']
        print(min(score),max(score))
        quant = df['FG.Quantity']
        log_quant = np.log10(quant)
        cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
        rt = df['EG.ApexRT']
        

    
    if 'DIANN' in name:
        df = df.dropna(subset=['Precursor.Quantity'])
        
        charge = df['Precursor.Charge']
        im = df['IM']
        localization = df['PTM.Site.Confidence']
        score = df['CScore']
        quant = df['Precursor.Quantity']
        log_quant = np.log10(quant)
        cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
        rt = df['RT']


    dict[name] = cleaned_log
    

# print(dict)
df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in dict.items() })
# print(df)

ax = df.plot.kde(color = )

ax

