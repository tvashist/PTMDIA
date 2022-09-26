import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


hybrid_library = pd.read_csv('Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro/20220803_134456_PTMDIAProject_DIACurveAnalysis_WithSpecLib_Report_addedFGLabel.tsv', delimiter= '\t',low_memory = False)
print('read')

SS_library = pd.read_csv('Z:/Helium_Tan/PTMDIAProject_PhosphoBGCurve/Outputs/SpectralLibSearch/Pro_SmallLibSearch_LocFilter/20220902_105134_PTMDIAProject_TimsTOFPro_DIACurveAnalysis_SmallLib0.75Loc_Report.tsv', delimiter= '\t',low_memory = False)
print('read')

heavies = pd.read_csv('Y:/LabMembers/Tan/DIA_QuantitativePTMs/Peptide_Lists/Modified_MatchesSpectronaut/Modified_Heavies.tsv', delimiter= '\t')['Modified_HeaviesAnnotated'][0:234]
print('read')



hybrid_heavies = hybrid_library.loc[hybrid_library['FG.LabeledSequence'].str.contains('Label')]        #Only look at the entries that contain heavy modification
SS_heavies = SS_library.loc[SS_library['FG.LabeledSequence'].str.contains('Label')]



#Change annotations to be integer annotations of mods including heavy labeling
def modify(df, column):
    new_annotations = []

    for seq in df[column]:  # This column contains sequences with all mod annotations including heavy labeling
        new = seq.replace('Phospho (STY)', '+80').replace('Carbamidomethyl (C)', '+57').replace('Label:13C(6)15N(4)', '+10').replace('Label:13C(6)15N(2)', '+8').replace('[Oxidation (M)]', '').replace('[Label:13C(5)15N(1)]', '').replace('[Label:13C(6)15N(2)]', '').replace('[Label:T-5]', '').replace('[Label:13C(6)15N(1)]', '').replace('[Acetyl (Protein N-term)]', '')
        new = new[1:-1]
        new_annotations.append(new)

    df['IntAnnotations'] = new_annotations

modify(hybrid_heavies, 'FG.LabeledSequence')
modify(SS_heavies, 'FG.LabeledSequence')

found_in_hybrid = []

print('Completed modifying both data frames, moving on to iterations')


for seq in heavies:
    if seq in set(hybrid_heavies['IntAnnotations']):
        found_in_hybrid.append(seq)


print('Starting to build dictionary of scores...')

scores = {}
for h in found_in_hybrid:

    if h not in scores:
        scores[h] = {}
        scores[h]['Hybrid Score'] = None
        scores[h]['Single Shot Score'] = None


    single_hybrid = hybrid_heavies.loc[hybrid_heavies['IntAnnotations'] == h]           #This is not the problem
    hybrid_score = single_hybrid['EG.Cscore'].mean()
    scores[h]['Hybrid Score'] = hybrid_score

    single_SS = SS_heavies.loc[SS_heavies['IntAnnotations'] == h]
    SS_score = single_SS['EG.Cscore'].mean()
    scores[h]['Single Shot Score'] = SS_score


print('Finished building dictionary')


scoreDF = pd.DataFrame.from_dict(scores)
scoreDF.to_csv("Y:/LabMembers/Tan/DIA_QuantitativePTMs/Analyses/ScoreComparison/Scores.tsv", sep = '\t')


















