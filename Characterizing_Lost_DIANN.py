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


#For DIANN
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['DIANN_Mods']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['DIANN_Mods']


# In[4]:


DIANN_DDALib = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\DDA_Library\Three_varmods\dia_nn\out\\report.tsv', sep = '\t')


# In[ ]:





# In[5]:


spiked = DIANN_DDALib.loc[DIANN_DDALib['Modified.Sequence'].isin(heavies)|DIANN_DDALib['Modified.Sequence'].isin(lights)]
highest = spiked.loc[spiked['Run'].str.contains('10p0')]
lowest = spiked.loc[spiked['Run'].str.contains('0p1')]

#What's lost along the way?
ten = highest['Modified.Sequence']
one = lowest['Modified.Sequence']

unique = set(ten) - set(one)

highest_only = highest.loc[highest['Modified.Sequence'].isin(unique)]
highest_shared = highest.loc[~highest['Modified.Sequence'].isin(unique)]


print(len(highest_only), len(highest_shared))


# In[19]:


for index, row in highest_only.iterrows():
    seq = row['Modified.Sequence']
    rt = row['RT']
    
    if rt < 6:
        print(seq)


# In[14]:


list = [highest_shared, highest_only]
names = ['Kept', 'Lost']

dict = {}

for i in range(0, len(list)):
    df = list[i]
    name = names[i]
    
    df = df.dropna(subset=['Precursor.Quantity'])
    
    sequence = df['Modified.Sequence']
    charge = df['Precursor.Charge']
    im = df['IM']
    localization = df['PTM.Site.Confidence']
    score = df['CScore']
    quant = df['Precursor.Quantity']
    log_quant = np.log10(quant)
    cleaned_log = [v for v in log_quant if not math.isnan(v) and not math.isinf(v)]
    rt = df['RT']
#     print(min(rt), max(rt))
    
    dict[name] = rt
    

df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in dict.items() })
for x in df["Lost"]:
    if x > 0:
        print(x)

ax = df.plot.kde()

ax

