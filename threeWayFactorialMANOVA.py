# -*- coding: utf-8 -*-
"""
Created on Mon May  3 12:41:00 2021

@author: Brian Chung
This script carries out three-way factorial MANOVAs on enzyme parameters across
3 independent variables: vegetation, precipitation, and time points (which
signify season). It will produce spreadsheets of MANOVA output (probably only
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
from statsmodels.multivariate.manova import MANOVA

# Loading in the log base 10 transformed Michaelis-Menten parameters to perform
# factorial ANOVA on
cwd = Path(os.getcwd())
cwdDirsFiles = os.listdir(cwd)
dataFolder = cwd/'Enzyme activity data'
dataDirsFiles = os.listdir(dataFolder)
paramsPath = dataFolder/'Parameters - log 10 transformed.xlsx'
parameters = pd.read_excel(paramsPath)
# %%
# Purpose: Performing a factorial MANOVA with the factors being time point,
# vegetation, and precipitation. This factorial MANOVA is used as exploratory
# data analysis to see which relationships are potentially significant.
'''Ok so just talked with Steve and he recommended the following tests in the
following order: factorial MANOVA for all enzymes, followed by factorial
ANOVAs & Tukey for individual enzymes. So Imma do that instead'''

# (1) Convert the parameters dataframe into a MultiIndex
parameters = parameters.drop(columns=["Transformation"])
indices = ["timePoint", "ID", "Vegetation", "Precip", "Replicate", "Enzyme",
           "Parameter"]
paramsMulInd = parameters.set_index(keys=indices)
paramsUnstack1 = paramsMulInd.unstack("Parameter")
paramsUnstack2 = paramsUnstack1.unstack("Enzyme")
indicesPivot = ["timePoint", "ID", "Vegetation", "Precip", "Replicate"]
indicesCols = ["Enzyme", "Parameter"]
paramsPivot = parameters.pivot(index=indicesPivot, columns=indicesCols,
                               values="value")
'''So I accomplished the same thing with both pivot and unstack, just that
they have slightly different levels with unstack resulting in the original
value column being an extra level compared to pivot'''

# (2) Making a wide version of the parameters dataframe to better suite
# the way statsmodels run MANOVAs
'''So my initial attempt at making a wide version of the parameters dataframe
resulted in NaNs all over the place that statsmodels do not like and so I
couldn't run my MANOVA. So I have to create 2 separate dataframes, 1 for
hydrolase, and 1 for oxidase, and then merge them together.
'''
params1Level = paramsPivot.copy()
params1Level.columns = params1Level.columns.map("_".join).str.strip("_")
params1Level.reset_index(inplace=True)
paramsHydro = params1Level.drop(columns=["Replicate", "PPO_Vmax", "PPO_Km"])
paramsHydro = paramsHydro.dropna()
hydroColsToDrop = ["AG_Vmax", "AP_Vmax", "BG_Vmax", "BX_Vmax", "CBH_Vmax",
                   "LAP_Vmax", "NAG_Vmax", "AG_Km", "AP_Km", "BG_Km", "BX_Km",
                   "CBH_Km", "LAP_Km", "NAG_Km"]

paramsOxi = params1Level.drop(columns=hydroColsToDrop)
paramsOxi = paramsOxi.dropna()
oxiIndices = ["timePoint", "ID", "Vegetation", "Precip"]
oxiPivotVals = ["PPO_Vmax", "PPO_Km"]
paramsOxi["Replicate"] = paramsOxi["Replicate"].astype(str)
paramsOxiPivot = paramsOxi.pivot(oxiIndices, "Replicate", oxiPivotVals)
paramsOxiPivot.columns = paramsOxiPivot.columns.map("_".join).str.strip("_")
paramsOxiPivot.columns = paramsOxiPivot.columns.str.strip(".0")
paramsOxiPivot = paramsOxiPivot.dropna(axis="columns")
paramsOxiPivot.reset_index(inplace=True)
paramsWide = pd.merge(paramsHydro, paramsOxiPivot, on=["timePoint", "ID",
                                                       "Vegetation", "Precip"])

# (3) Performing the MANOVA analysis
manovaFormula = "AG_Vmax + AG_Km + AP_Vmax + AP_Km + BG_Vmax + BG_Km + BX_Vmax + BX_Km + CBH_Vmax + CBH_Km + LAP_Vmax + LAP_Km + NAG_Vmax + NAG_Km + PPO_Vmax_1 + PPO_Vmax_2 + PPO_Km_1 + PPO_Km_2 ~ C(timePoint, Sum)*C(Vegetation, Sum)*C(Precip, Sum)"
manovaModel = MANOVA.from_formula(formula=manovaFormula, data=paramsWide)
manovaResults = manovaModel.mv_test()
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
initFormula = "value ~ C(timePoint, Sum)*C(Vegetation, Sum)*C(Precip, Sum)"
initModel = ols(formula=initFormula, data=VmaxAG).fit()
initialResults = anova_lm(initModel, typ=3)
'''So the three-way interaction isn't significant, time to split it up into
the 2-way interactions'''

# (3) Splitting up three-way interaction
timeXveg = "C(timePoint, Sum)*C(Vegetation, Sum)"
timeXppt = "C(timePoint, Sum)*C(Precip, Sum)"
vegXppt = "C(Vegetation, Sum)*C(Precip, Sum)"
sepInteracForm = "value ~ {0} + {1} + {2}".format(timeXveg, timeXppt, vegXppt)
sepInteracModel = ols(formula=sepInteracForm, data=VmaxAG).fit()
results2 = anova_lm(sepInteracModel, typ=3)
'''So there are no interaction effects between any of the 3 independent
variables on AG Vmax parameter values, even when they're all separated into
2-way interactions. So, let's drop all the non-significant interactions'''

# (4) ANOVA with no interactions
noInteracForm = "value ~ C(timePoint) + Vegetation + Precip"
noInteracModel = ols(formula=noInteracForm, data=VmaxAG).fit()
results3 = anova_lm(noInteracModel, typ=2)
'''So the main effects alone (time point, vegetation, and precip) don't have
any effect on AG either.'''
