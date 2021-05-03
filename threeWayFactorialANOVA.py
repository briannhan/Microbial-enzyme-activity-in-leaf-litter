# -*- coding: utf-8 -*-
"""
Created on Mon May  3 12:41:00 2021

@author: Brian Chung
This script carries out three-way factorial ANOVAs on enzyme parameters across
3 independent variables: vegetation, precipitation, and time points (which
signify season). It will produce spreadsheets of ANOVA output (probably only
1 spreadsheet for combinations of vegetation, precipitation, and time points,
and 1 spreadsheet for combinations of vegetation, precipitation, and season,
if I choose to classify time points by season).
"""
import pandas as pd
import os
from pathlib import Path
# import statsmodels as sm
# import statsmodels.formula.api as smfa
from statsmodels.formula.api import ols
# from statsmodels.stats.multicomp import MultiComparison
from statsmodels.stats.anova import anova_lm

# Loading in the log base 10 transformed Michaelis-Menten parameters to perform
# factorial ANOVA on
cwd = Path(os.getcwd())
cwdDirsFiles = os.listdir(cwd)
dataFolder = cwd/'Enzyme activity data'
dataDirsFiles = os.listdir(dataFolder)
paramsPath = dataFolder/'Parameters - log 10 transformed.xlsx'
parameters = pd.read_excel(paramsPath)
# %%
# Purpose: Carry out three-way ANOVA on AG's Vmax as a proof of concept to code
# ANOVAs for the remaining 7 enzymes. This process will be run iteratively if
# interactions between at least 2 independent variables (vegetation, time
# points, and precipitation) are non-significant, and each non-significant
# interaction will be broken down into simpler interactions and removed until
# the significant simpler interactions remain.

# (1) Subsetting out AG's Vmax
conditions = (parameters.Enzyme == "AG") & (parameters.Parameter == "Vmax")
VmaxAG = parameters[conditions]

# (2) Initial linear regression model and ANOVA attempt
initFormula = "value ~ C(timePoint, Sum) + C(Vegetation, Sum) + C(Precip, Sum) + C(timePoint, Sum)*C(Vegetation, Sum)*C(Precip, Sum)"
initModel = ols(formula=initFormula, data=VmaxAG).fit()
initialResults = anova_lm(initModel, typ=3)
'''So the three-way interaction isn't significant, time to split it up into
the 2-way interactions'''

# (3) Splitting up three-way interaction
sepInteractionsForm = "value ~ C(timePoint, Sum) + C(Vegetation, Sum) + C(Precip, Sum) + C(timePoint, Sum):C(Vegetation, Sum) + C(timePoint, Sum):C(Precip, Sum) + C(Vegetation, Sum):C(Precip, Sum)"
sepInteractionsModel = ols(formula=sepInteractionsForm, data=VmaxAG).fit()
results2 = anova_lm(sepInteractionsModel, typ=3)
'''So there are no interaction effects between any of the 3 independent
variables on AG Vmax parameter values, even when they're all separated into
2-way interactions. So, let's drop all the non-significant interactions'''

# (4) ANOVA with no interactions
noInteractionsForm = "value ~ C(timePoint, Sum) + C(Vegetation, Sum) + C(Precip, Sum)"
noInteractionsModel = ols(formula=noInteractionsForm, data=VmaxAG).fit()
results3 = anova_lm(noInteractionsModel, typ=3)
'''So the main effects alone (time point, vegetation, and precip) don't have
any effect on AG either.'''
