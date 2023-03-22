#!/usr/bin/env python
# coding: utf-8

# In[1]:


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[16]:


import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


# In[20]:


diann = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\DDA_Library\Three_varmods\dia_nn\out\\report.tsv', sep = '\t', low_memory=False)
# diann = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Combined_DDA_Predicted_Library\dia_nn\out\\report.tsv', sep = '\t')
# diann = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\DIANN\Predicted_Library\dia_nn\out\\report.tsv', sep = '\t')

print('read')


# In[21]:


platform = 'Pro'

spikes = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
spikes = [str(x).replace('.','p') for x in spikes]

if platform == 'SCP':
    conditions = spikes[:9]
else:
    conditions = spikes[4:]
    
print(conditions)


# In[30]:


#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['DIANN_Mods']
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter = '\t')['DIANN_Mods'][0:234]


singly_phosphorylated = [x for x in list(heavies) if str(x).count('(UniMod:21)')==1]
# singly_phosphorylated.extend(list(lights)) 

print(len(heavies),len(singly_phosphorylated) )


# In[1]:


# def scoring():
    
#     scoring = {}
    
    
#     for c in conditions:
#         conc = diann.loc[diann['Run'].str.contains(str(c))]
#         scores = probabilities = list(conc['CScore'])
        
#         scoring[c] = scores
        
#     df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in scoring.items() })
#     print(df)
# #     ax = df.plot.kde()
#     df['0p1'].hist()
#     df['10p0'].hist()
    

# scoring()
    
    


# In[31]:


def localize():

    dist = {}
    localized = {}

    for c in conditions:

        probabilities = []

        
        conc = diann.loc[diann['Run'].str.contains(str(c))]                                #Only runs from that spiked condition
        spiked = conc.loc[conc['Modified.Sequence'].isin(singly_phosphorylated)]           #Only spiked peptide entries
        probabilities = list(spiked['PTM.Site.Confidence'])

        a = c.replace('p','.')
        dist[a] = probabilities
            
        
        reg = [x for x in probabilities if x >= 0.01]
        stringent = [x for x in probabilities if x >= 0.51]

        
        perc_reg = (len(reg)/len(probabilities))*100
        perc_stringent = (len(stringent)/len(probabilities))*100
        
        
        
        localized[a] = (perc_reg, perc_stringent)
        print(c, len(probabilities),perc_reg, perc_stringent)
            
    

    df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in dist.items() })
    localized_df = pd.DataFrame(localized, index = ['0.01','0.51'])
    localized_df = localized_df.transpose()



    ax = df.plot.kde()
    bx = localized_df.plot.bar()







localize()


# In[40]:


dict = {}
reg = [x for x in range(0,100)]
stringent = [x for x in range(0,100)]*6

dict['r'] = reg
dict['s'] = stringent

df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in dict.items() })

print(df)
ax = df.plot.kde()

