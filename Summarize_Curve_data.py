import pandas as pd
import re
import math
from statistics import mean,stdev
import os

platform = 'Exploris'

if platform == 'SCP':
    conditions = ['0.004fmol', '0.01fmol', '0.02fmol', '0.04fmol', '0.1fmol', '0.2fmol', '0.4fmol', '1.0fmol', '2.0fmol', '4.0fmol']


if platform == 'Pro' or platform == 'Exploris':
    conditions = ['0.04fmol', '0.1fmol', '0.2fmol', '0.4fmol', '1.0fmol', '2.0fmol', '4.0fmol', '10fmol']


# assign directory
directory = 'Exploris_Lights_outputs_Found'

# iterate over files in that directory


collect_reps = []
collect_cv = []

for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    if 'tsv' in f:                          #Only tsv files, not pdf images
        df = pd.read_csv(f, delimiter= '\t')



