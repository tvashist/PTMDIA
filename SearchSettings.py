import pandas as pd
import numpy as np
import re
import math
from statistics import mean,stdev,mode

#First: Convert modifications as annotated in Spectronaut results ouput to simpler annotations including for heavies

setting = '5TransNoLoc'
spectronaut = pd.read_csv('S:/Helium_Tan/PTMDIAProject_SearchForHeavies/Out/WithPhosphoBG/5TransNoLoc/20220728_142424_PTMDIAProject_SpecLibTest_5TransNoLocFilter_HeaviesPhosphoBG_Report.tsv', delimiter= '\t')
read_heavies = pd.read_csv('Z:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Heavies.tsv', delimiter= '\t')
heavies = list(read_heavies['Modified_HeaviesAnnotated'])
heavies = [x for x in heavies if pd.isnull(x) == False]


new_annotations = []
for seq in spectronaut['FG.LabeledSequence']:      #This column contains sequences with all mod annotations including heavy labeling
    new = seq.replace('Phospho (STY)', '+80').replace('Carbamidomethyl (C)','+57').replace('Label:13C(6)15N(4)','+10').replace('Label:13C(6)15N(2)','+8').replace('[Oxidation (M)]', '+16')
    new = new[1:-1]
    new_annotations.append(new)
spectronaut['AnnotatedSequence'] = new_annotations
sub = spectronaut[['R.Condition','R.Replicate','EG.ModifiedSequence','FG.Charge', 'FG.LabeledSequence','AnnotatedSequence']]

count = 0
for h in heavies:
    # print(h)
    if h in sub['AnnotatedSequence'].values:
        count += 1

print('The number of heavies detected searching against the library with settings ' + setting + ': ' + str(count))




