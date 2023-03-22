#!/usr/bin/env python
# coding: utf-8

# In[20]:


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[21]:


import pandas as pd
import matplotlib
from matplotlib_venn import venn2, venn3
from upsetplot import plot
import seaborn as sns
import numpy as np
import math


# In[8]:


SN_DDALib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\DDALibrary\\20230308_091202_R02_PhosphoDIA_Pro_DDALibrary_Report.tsv', sep = '\t')


# In[23]:


FilteredOut = SN_DDALib.loc[SN_DDALib['EG.PTMAssayProbability'] < 0.75]
s = FilteredOut['EG.PTMAssayProbability']


ax = s.plot.kde()
ax


# In[9]:


lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['LibraryHeavies'][0:234]


# In[22]:


print(SN_DDALib['R.Condition'].unique())


# In[14]:


spiked = SN_DDALib.loc[SN_DDALib['FG.IntMID'].isin(heavies)|SN_DDALib['FG.IntMID'].isin(lights)]
highest = spiked.loc[spiked['R.Condition'].str.contains('10.0')]
lowest = spiked.loc[spiked['R.Condition'].str.contains('0.1')]

#What's lost along the way?
ten = highest['FG.IntMID']
one = lowest['FG.IntMID']

unique = set(ten) - set(one)

highest_only = highest.loc[highest['FG.IntMID'].isin(unique)]
highest_shared = highest.loc[~highest['FG.IntMID'].isin(unique)]


print(len(highest_only), len(highest_shared))


# In[24]:


list = [highest_shared, highest_only]
names = ['Kept', 'Lost']

dict = {}

for i in range(0, len(list)):
    df = list[i]
    name = names[i]
    
    df = df.dropna(subset=['EG.PTMAssayProbability'])


    charge = df['FG.Charge']
    im = df['FG.IonMobility']
    localization =df['EG.PTMAssayProbability']
    score = df['EG.Cscore']
    quant = df['FG.Quantity']
    log_quant = np.log10(quant)
    cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
    rt = df['EG.ApexRT']
    
    dict[name] = rt
    


df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in dict.items() })

ax = df.plot.kde()

ax

