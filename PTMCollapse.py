import pandas as pd
import re
import math
from statistics import mean,stdev,mode
import os


wd = 'Y:\\LabMembers\\Tan\\CompGroup\\PTMSiteCollapse\\'
file = '20221122_163329_PTMDIAProject_ExplorisFAIMS_3SSLib_DIACurveAnalysis_Report.tsv'

report = pd.read_csv(wd + file, delimiter ='\t', low_memory= False)

short = report[0:100]
short.to_csv(wd + 'short.tsv', sep = '\t', index = False)

