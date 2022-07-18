import pandas as pd

lights = pd.read_csv('NHPeptideList.csv')

modified = []
for seq in lights['Peptides']:
    for x in seq:
        if x.islower():
            upper = x.upper()
            new = seq.replace(x, upper+'[+80]')
            final = '_'+new+'_'
            modified.append(final)

lights['Modified'] = modified
lights.to_csv('Modified_Lights.tsv', sep = '\t')


heavies = pd.read_csv('HuHeavypSTY_PeptideList.csv')

modified_heavy = []
for seq in heavies['Phosphopeptide Sequence']:
    new = seq.replace('(Cam)', '').replace('(pS)', 'S[+80]').replace('(pT)', 'T[+80]').replace('(pY)', 'Y[+80]')
    final = '_'+new+'_'
    modified_heavy.append(final)

heavies['Modified'] = modified_heavy
heavies = heavies[['Phosphopeptide Sequence', 'Modified']]

heavies.to_csv('Modified_Heavies.tsv', sep = '\t')








