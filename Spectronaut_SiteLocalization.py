#!/usr/bin/env python
# coding: utf-8

# In[2]:


from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[3]:


import pandas as pd
import re


# In[12]:


# spectronaut = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\directDIA\\20230227_104026_PhosphoDIA_R02_Report_directDIA.tsv', sep = '\t', low_memory = False)
spectronaut = pd.read_csv('Z:\Helium_Tan\R02_PTMDIA\Pro\Spectronaut\Hybrid\\20230313_103651_R02_PhosphoDIA_Pro_HybridLibrary_Report.tsv', sep = '\t', low_memory = False)
print('read')


# In[13]:


platform = 'Pro'

spikes = [0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0]
if platform == 'SCP':
    conditions = spikes[:9]
else:
    conditions = spikes[4:]
    
print(conditions)


# In[10]:


#Modified lights document in lab members folder
lights = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Lights.tsv', delimiter= '\t')['Modified']        #Modified sequence annotates phosphopeptides as they appear on Spectronaut
heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified/Modified_Heavies.tsv', delimiter= '\t')['LibraryHeavies']


singly_phosphorylated = [x for x in list(heavies) if str(x).count('[+80]')==1]   #Only singly phosphorylated heavies
singly_phosphorylated.extend(list(lights))                                       #Add non-human phosphopeptides, which are all singly-phosphorylated


# In[14]:


def localize():

    dist = {}
    localized = {}


    for c in conditions:

        probabilities = []

        conc = spectronaut.loc[spectronaut['R.Condition'].str.contains(str(c))]                                #Only runs from that spiked condition
        spiked = conc.loc[conc['FG.IntMID'].isin(singly_phosphorylated)]
        print(len(spiked))

        for index, row in spiked.iterrows():
            
            probs = row['EG.PTMProbabilities [Phospho (STY)]']
            
            if probs != float:
                prob = probs.split(';')

            seq = row['PEP.StrippedSequence']
            mod_seq = row['EG.IntPIMID']

            #Finding probabilitiy info

            all_sty = []
            find_sty = re.finditer('[STY](\[\+80\])?', mod_seq)    #Find position of S|T|Y, regardless of if followed by +80 or not
            for match in find_sty:                       
                all_sty.append(match.start())

            all_phos = []
            find_phos = re.finditer('[STY](\[\+80\])', mod_seq)    #Find position of phosphorylated S|T|Y, which are followed by +80
            for match in find_phos:
                all_phos.append(match.start())

            index = [all_sty.index(x) for x in all_phos][0]
            actual_prob = prob[index]
            
            
            probabilities.append(float(actual_prob))
        
            
        
        reg = [x for x in probabilities if x >= 0.75]
        stringent = [x for x in probabilities if x >= 0.99]
        
        perc_reg = (len(reg)/len(probabilities))*100
        perc_stringent = (len(stringent)/len(probabilities))*100
        
        
        
        localized[c] = (perc_reg, perc_stringent)
            
        
        dist[c] = probabilities
        
    df = pd.DataFrame.from_dict({ key:pd.Series(value) for key, value in dist.items() })
    localized_df = pd.DataFrame(localized, index = ['0.75','0.99'])
    localized_df = localized_df.transpose()



    ax = df.plot.kde()
    bx = localized_df.plot.bar()







localize()

