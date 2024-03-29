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

Analysis is heavily based on the following links:
https://www.pythonfordatascience.org/factorial-anova-python/
https://youtu.be/d_Azlncd-kU
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
from datetime import datetime
start = datetime.now()

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
"""The Sum contrast allows for type 1, 2, and 3 ANOVAs to have the exact same
sum of squares for the main effects, their interactions, and the residuals. If
Sum is not included in the formula, then type 1 & 2 sum of squares will still
be the same as if Sum was used, but type 3 sum of squares won't be the same. In
other words, only type 1 & 2 sum of squares will be correctly calculated."""
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


# (5) Writing a function to abstract the process of performing a factorial
# ANOVA
def factorialANOVA(enzyme, parameter, threeWay=None, twoWay1=None, twoWay2=None, twoWay3=None):
    """
    Carries out a three-way ANOVA on a specific parameter of an enzyme. Each
    time it is called, it runs an ANOVA once.

    Parameters
    ----------
    enzyme : str
        The name of the enzyme whose parameters will undergo a factorial ANOVA.
    parameter : str
        The enzyme parameter that will undergo a factorial ANOVA.
    threeWay : str, optional
        The three-way interaction, if specified, and the user wants to test for
        a three-way interaction. If specified, then a type 3 ANOVA will be
        performed. The value of this argument must be 'Y' for 'yes' if the user
        wants to test for this type of interaction. Otherwise, if the user does
        not want to test for a three-way interaction due to a previous test
        revealing that there's no significant three-way interaction, then this
        can be left blank. The default is None.
    twoWay1 : str, optional
        A two-way interaction that the user wants to test for. If this
        argument, along with all the other interactions arguments, is left
        blank, then this means that there is no two-way interactions at all
        and a type 2 ANOVA will be performed. The default is None.
    twoWay2 : str, optional
        A two-way interaction that the user wants to test for. If this
        argument, along with all the other interactions arguments, is left
        blank, then this means that there is no two-way interactions at all
        and a type 2 ANOVA will be performed. The default is None.
    twoWay3 : str, optional
        A two-way interaction that the user wants to test for. If this
        argument, along with all the other interactions arguments, is left
        blank, then this means that there is no two-way interactions at all
        and a type 2 ANOVA will be performed. The default is None.

    Returns
    -------
    Pandas dataframe of ANOVA results

    """
    boolean1 = parameters["Enzyme"] == enzyme
    boolean2 = parameters["Parameter"] == parameter
    dataParams = parameters[(boolean1) & (boolean2)]
    if threeWay == "Y":  # Testing for three-way interactions
        formu = "value ~ C(timePoint, Sum)*C(Vegetation, Sum)*C(Precip, Sum)"
        model = ols(formula=formu, data=dataParams).fit()
        results = anova_lm(model, typ=3)
    elif twoWay1 is not None:  # Testing for at least 1 two-way interaction
        if twoWay2 is not None:  # Testing for at least 2 two-way interactions
            if twoWay3 is not None:  # Testing for 3 two-way interactions
                formu = "value ~ {0} + {1} + {2}".format(twoWay1, twoWay2,
                                                         twoWay3)
            else:
                formu = "value ~ {0} + {1}".format(twoWay1, twoWay2)
        else:
            factor1 = "C(timePoint, Sum)"
            factor2 = "C(Vegetation, Sum)"
            factor3 = "C(Precip, Sum)"
            formu = "value ~ {0} + {1} + {2} + {3}".format(factor1, factor2,
                                                           factor3, twoWay1)
        model = ols(formula=formu, data=dataParams).fit()
        results = anova_lm(model, typ=3)
    else:  # If there are no interactions, so a Type 2 ANOVA will be ran
        formu = "value ~ C(timePoint) + Vegetation + Precip"
        model = ols(formula=formu, data=dataParams).fit()
        # if enzyme != "PPO":
        #     results = anova_lm(model, typ=2)
        # elif enzyme == "PPO":
        #     # PPO undergoes a type 3 ANOVA
        #     results = anova_lm(model, typ=3)
        results = anova_lm(model, typ=2)

    # Calculating eta-squared and partial eta-squared
    residualSS = results.loc["Residual", "sum_sq"]
    for index, row in results.iterrows():
        rowSS = row.sum_sq
        results.loc[index, "partialEta2"] = rowSS/(residualSS + rowSS)
    return results


# %%
# Purpose: Running factorial ANOVAs for AG Km

# (1) Testing for all interactions (two-way and three-way)
results1KmAG = factorialANOVA(enzyme="AG", parameter="Km", threeWay="Y")
'''Three-way interactions are not significant, so removing that and specifying
the individual two-way interactions'''

# (2) Testing for two-way interactions only
results2KmAG = factorialANOVA(enzyme="AG", parameter="Km", twoWay1=timeXveg,
                              twoWay2=timeXppt, twoWay3=vegXppt)
'''No two-way interactions are significant either, so now I'm gonna test for
only the main effects.'''

# (3) Testing for main effects only
results3KmAG = factorialANOVA(enzyme="AG", parameter="Km")
'''And there are no significant main effects. So looks like time point,
vegetation, and precipitation does not influence alpha-glucosidase activity,
meaning that alpha-glucosidase activity is relatively constant between both
coastal sage scrub and grasslands regardless of precipitation treatments, the
amount of time that decomposition had occurred, and season (dry or growing).'''
# %%
# Purpose: Running factorial ANOVAs for AP Vmax

# (1) Testing for all interactions (two-way and three-way)
results1VmaxAP = factorialANOVA("AP", "Vmax", "Y")
'''So timePoint x vegetation and timePoint x precipitation are almost
significant, but otherwise no significant interactions. timePoint as a main
effect is signifcant though. Time to remove the three-way interaction.'''

# (2) Testing for two-way interactions only
results2VmaxAP = factorialANOVA("AP", "Vmax", twoWay1=timeXppt,
                                twoWay2=timeXveg, twoWay3=vegXppt)
'''Again, timePoint x vegetation and timePoint x precipitation are only
marginally significant, while the other interaction,
vegetation x precipitation, isn't significant. Also, again, timePoint as a main
effect is significant, but this time, precipitation as a main effect is
marginally significant. Time to remove all the two-way interactions.'''

# (3) Testing for main effects only
results3VmaxAP = factorialANOVA("AP", "Vmax")
'''timePoint as a main effect is significant, while precipitation as a main
effect is only marginally significant this time around.'''

# %%
# Purpose: Running factorial ANOVAs for AP Km

# (1) Testing for all interactions (two-way and three-way)
results1KmAP = factorialANOVA("AP", "Km", "Y")

# (2) Testing for two-way interactions only
results2KmAP = factorialANOVA("AP", "Km", twoWay1=timeXppt, twoWay2=timeXveg,
                              twoWay3=vegXppt)

# (3) Removing timePoint x precipitation, which is the only non-signifcant
# two-way interaction
results3KmAP = factorialANOVA("AP", "Km", twoWay1=timeXveg, twoWay2=vegXppt)
# %%
# Purpose: Running factorial ANOVAs for BG Vmax

# (1) Testing all interactions (two-way and three-way)
results1VmaxBG = factorialANOVA("BG", "Vmax", "Y")

# (2) Testing only two-way interactions
results2VmaxBG = factorialANOVA("BG", "Vmax", twoWay1=timeXppt,
                                twoWay2=timeXveg, twoWay3=vegXppt)

# (3) Testing model with no interactions
results3VmaxBG = factorialANOVA("BG", "Vmax")

# (4) Testing whether I should be wrapping vegetation and precip in C()
BGparams = parameters[(parameters.Enzyme == "BG")
                      & (parameters.Parameter == "Vmax")]
formuVegPptWrap = "value ~ C(timePoint) + C(Vegetation) + C(Precip)"
modelVegPptWrap = ols(formula=formuVegPptWrap, data=BGparams).fit()
resultsVegPptWrap = anova_lm(modelVegPptWrap, typ=2)
"""So it makes no difference whether or not I wrap Vegetation or Precip in C()
"""

# (5) Testing what happens if I don't wrap timePoint in C()
formuTimeUnwrapped = "value ~ timePoint + Vegetation + Precip"
modelTimeUnwrapped = ols(formula=formuTimeUnwrapped, data=BGparams).fit()
resultsTimeUnwrapped = anova_lm(modelTimeUnwrapped, typ=2)
"""Not wrapping time, makes quite a big difference. I'm assuming that this is
because time is already a numerical variable (integer), so by not wrapping it,
the linear regression model just performs a regression on time instead of
treating time as categorical. The C() is used to explicitly denote a variable
as categorical. Since Vegetation and Precip have values of string in the first
place, they're categorical by default, but time isn't a string by default, so
I have to wrap time with C()."""
# %%
# Purpose: Running factorial ANOVAs for BG Km

# (1) Testing all interactions (two-way and three-way)
results1KmBG = factorialANOVA("BG", "Km", "Y")

# (2) Testing only two-way interactions
results2KmBG = factorialANOVA("BG", "Km", twoWay1=timeXppt, twoWay2=timeXveg,
                              twoWay3=vegXppt)

# (3) Testing main effects only
results3KmBG = factorialANOVA("BG", "Km")
# %%
# Purpose: Running factorial ANOVAs for BX Vmax

# (1) Testing all interactions (two-way and three-way)
results1VmaxBX = factorialANOVA("BX", "Vmax", "Y")

# (2) Testing only two-way interactions
results2VmaxBX = factorialANOVA("BX", "Vmax", twoWay1=timeXppt,
                                twoWay2=timeXveg, twoWay3=vegXppt)

# (3) Testing main effects only
results3VmaxBX = factorialANOVA("BX", "Vmax")

# %%
# Purpose: Running factorial ANOVAs for BX Km

# (1) Testing for all interactions (two-way and three-way)
results1KmBX = factorialANOVA("BX", "Km", "Y")

# (2) Testing two-way interactions only
results2KmBX = factorialANOVA("BX", "Km", None, timeXppt, timeXveg, vegXppt)

# (3) Testing main effects only
results3KmBX = factorialANOVA("BX", "Km")

# %%
# Purpose: Running factorial ANOVAs for CBH Vmax

# (1) Testing for all interactions (two-way and three-way)
results1VmaxCBH = factorialANOVA("CBH", "Vmax", "Y")

# (2) Testing for two-way interactions only
results2VmaxCBH = factorialANOVA("CBH", "Vmax", twoWay1=timeXppt,
                                 twoWay2=timeXveg, twoWay3=vegXppt)

# (3) Testing after removal of vegetation x precipitation
results3VmaxCBH = factorialANOVA("CBH", "Vmax", None, timeXppt, timeXveg)
# %%
# Purpose: Running factorial ANOVAs for CBH Km

# (1) Testing for all interactions (two-way and three-way)
results1KmCBH = factorialANOVA("CBH", "Km", "Y")

# (2) Testing two-way interactions only
results2KmCBH = factorialANOVA("CBH", "Km", None, timeXppt, timeXveg, vegXppt)

# (3) Removing vegetation x precipitation and testing for the remaining
# two-way interactions
results3KmCBH = factorialANOVA("CBH", "Km", None, timeXppt, timeXveg)
# %%
# Purpose: Running factorial ANOVAs for LAP Vmax

# (1) Testing for all interactions (two-way and three-way)
results1VmaxLAP = factorialANOVA("LAP", "Vmax", "Y")

# (2) Testing for two-way interactions only
results2VmaxLAP = factorialANOVA("LAP", "Vmax", twoWay1=timeXppt,
                                 twoWay2=timeXveg, twoWay3=vegXppt)

# (3) Testing for main effects only
results3VmaxLAP = factorialANOVA("LAP", "Vmax")
# %%
# Purpose: Running factorial ANOVAs for LAP Km

# (1) Testing for all interactions (two-way and three-way)
results1KmLAP = factorialANOVA("LAP", "Km", "Y")

# (2) Testing for two-way interactions only
results2KmLAP = factorialANOVA("LAP", "Km", None, timeXppt, timeXveg, vegXppt)

# (3) Testing for main effects only
results3KmLAP = factorialANOVA("LAP", "Km")
# %%
# Purpose: Running factorial ANOVAs for NAG Vmax

# (1) Testing for all interactions (two-way and three-way)
results1VmaxNAG = factorialANOVA("NAG", "Vmax", "Y")

# (2) Testing two-way interactions only
results2VmaxNAG = factorialANOVA("NAG", "Vmax", twoWay1=timeXppt,
                                 twoWay2=timeXveg, twoWay3=vegXppt)

# (3) Removing vegetation x precipitation and timePoint x precipitation, which
# are the non-significant two-way interactions, and testing for the remaining
# timePoint x vegetation interaction
results3VmaxNAG = factorialANOVA("NAG", "Vmax", None, timeXveg)
# %%
# Purpose: Running factorial ANOVAs for NAG Km

# (1) Testing for all interactions (two-way and three-way)
results1KmNAG = factorialANOVA("NAG", "Km", "Y")

# (2) Testing two-way interactions only
results2KmNAG = factorialANOVA("NAG", "Km", None, timeXppt, timeXveg, vegXppt)

# (3) Removing vegetation x precipitation and timePoint x precipitation, which
# are the non-significant two-way interactions, and testing for the remaining
# timePoint x vegetation interaction
results3KmNAG = factorialANOVA("NAG", "Km", None, timeXveg)
# %%
# Purpose: Running factorial ANOVAs for PPO Vmax

# (1) Testing for all interactions (two-way and three-way)
results1VmaxPPO = factorialANOVA("PPO", "Vmax", "Y")
'''Lol and I already got a significant three-way interaction right off the bat
'''
# %%
# Purpose: Running factorial ANOVAs for PPO Km

# (1) Testing for all interactions (two-way and three-way)
results1KmPPO = factorialANOVA("PPO", "Km", "Y")

# (2) Testing two-way interactions only
results2KmPPO = factorialANOVA("PPO", "Km", None, timeXppt, timeXveg, vegXppt)

# (3) Removing non-significant two-way interactions, which are
# vegetation x precipitation and timePoint x precipitation, and testing for the
# remaining significant two-way interaction, timePoint x vegetation
results3KmPPO = factorialANOVA("PPO", "Km", None, timeXveg)
# %%
# Purpose: Loading in and wrangling the litter chemistry data to also perform
# ANOVAs

# (1) Loading in litter chemistry data
litterChemFolder = cwd/"Litter chemistry"
litterChemPath = litterChemFolder/"Carbohydrates and Proteins FTIR.xlsx"
litterChem = pd.read_excel(litterChemPath)

# (2) Wrangling the litter chemistry data to have the same columns as the
# enzyme parameters dataframe to allow the litter chemistry data to be tested
# using the same ANOVA code
litterChem.rename({"id": "ID", "functionalGroup": "Enzyme"},
                  axis="columns", inplace=True)
# "functionalGroup" is renamed as "Enzyme" so that this dataframe can be merged
# with the Michaelis-Menten parameters dataframe and ANOVAs can be performed
carbohydrate = ["glycosidicBond", "C_O_stretching", "carboEster",
                "carboEster1", "carboEster2"]
protein = ["amide", "amide1", "amide2"]
for index, row in litterChem.iterrows():
    functionalGroup = row.Enzyme
    if functionalGroup in carbohydrate:
        litterChem.loc[index, "Parameter"] = "Carbohydrate"
    elif functionalGroup in protein:
        litterChem.loc[index, "Parameter"] = "Protein"
    elif functionalGroup == "lipid" or functionalGroup == "alkane":
        litterChem.loc[index, "Parameter"] = functionalGroup

# (3) Merging the litter chemistry and Michaelis-Menten parameter dataframes
# into a single dataframe
parameters = pd.concat([parameters, litterChem], join="outer")
# %%
# Purpose: Running factorial ANOVAs for glycosidic bonds

# (1) Testing for all interactions (two-way and three-way)
results1glycosidic = factorialANOVA("glycosidicBond", "Carbohydrate", "Y")

# (2) Testing two-way interactions only
results2glycosidic = factorialANOVA("glycosidicBond", "Carbohydrate", None,
                                    timeXppt, timeXveg, vegXppt)

# (3) Removing two-way interactions with timePoint, leaving only the
# Vegetation x Precipitation interaction in
results3glycosidic = factorialANOVA("glycosidicBond", "Carbohydrate", None,
                                    vegXppt)
# %%
# Purpose: Running factorial ANOVAs for carbohydrate C_O_stretching

# (1) Testing for all interactions (two-way and three-way)
results1C_O = factorialANOVA("C_O_stretching", "Carbohydrate", "Y")

# (2) Testing for two-way interactions only
results2C_O = factorialANOVA("C_O_stretching", "Carbohydrate", None, timeXppt,
                             timeXveg, vegXppt)

# (3) Removing all two-way interactions, which were all insignificant
results3C_O = factorialANOVA("C_O_stretching", "Carbohydrate")
# %%
# Purpose: Running factorial ANOVAs for total carbohydrate esters, i.e.
# carboEster1 + carboEster2

# (1) Testing for all interactions (two-way and three-way)
results1carboEster = factorialANOVA("carboEster", "Carbohydrate", "Y")

# (2) Testing for two-way interactions only
results2carboEster = factorialANOVA("carboEster", "Carbohydrate", None,
                                    timeXppt, timeXveg, vegXppt)

# (3) Testing for main effects only, since all two-way interactions are
# insignificant
results3carboEster = factorialANOVA("carboEster", "Carbohydrate")
# %%
# Purpose: Running factorial ANOVAs for carboEster1

# (1) Testing for all interactions (two-way and three-way)
results1carboEster1 = factorialANOVA("carboEster1", "Carbohydrate", "Y")

# (2) Testing for two-way interactions only
results2carboEster1 = factorialANOVA("carboEster1", "Carbohydrate", None,
                                     timeXppt, timeXveg, vegXppt)

# (3) Testing for main effects only, since all two-way interactions are
# insignificant
results3carboEster1 = factorialANOVA("carboEster1", "Carbohydrate")
# %%
# Purpose: Running factorial ANOVAs for carboEster2

# (1) Testing for all interactions (two-way and three-way)
results1carboEster2 = factorialANOVA("carboEster2", "Carbohydrate", "Y")

# (2) Testing for two-way interactions only, since the three-way interaction is
# insignificant
results2carboEster2 = factorialANOVA("carboEster2", "Carbohydrate", None,
                                     timeXppt, timeXveg, vegXppt)

# (3) Testing for only the veg x precip two-way interaction plus main effects,
# since other two-way interactions are insignificant
results3carboEster2 = factorialANOVA("carboEster2", "Carbohydrate", None,
                                     vegXppt)
# %%
# Purpose: Running factorial ANOVAs for lipid content

# (1) Testing for all interactions (two-way and three-way)
results1lipid = factorialANOVA("lipid", "lipid", "Y")

# (2) Testing for two-way interactions only
results2lipid = factorialANOVA("lipid", "lipid", None, timeXppt, timeXveg,
                               vegXppt)

# (3) Testing for main effects and only the vegetation x precipitation two-way
# interaction
results3lipid = factorialANOVA("lipid", "lipid", None, vegXppt)
# %%
# Purpose: Running factorial ANOVAs on alkanes

# (1) Testing for all interactions (two-way and three-way)
results1alkane = factorialANOVA("alkane", "alkane", "Y")

# (2) Testing for two-way interactions only
results2alkane = factorialANOVA("alkane", "alkane", None, timeXppt, timeXveg,
                                vegXppt)

# (3) Testing for main effects and the two-way interactions of
# time x vegetation and precipitation x vegetation
results3alkane = factorialANOVA("alkane", "alkane", None, timeXveg, vegXppt)
# %%
# Purpose: Running factorial ANOVAs on total amides (amide1 + amide2)

# (1) Testing for all interactions (two-way and three-way)
results1protein = factorialANOVA("amide", "Protein", "Y")

# (2) Testing for two-way interactions only
results2protein = factorialANOVA("amide", "Protein", None, timeXppt, timeXveg,
                                 vegXppt)

# (3) Testing for main effects only
results3protein = factorialANOVA("amide", "Protein")
# %%
# Purpose: Running factorial ANOVAs on amide1

# (1) Testing for all interactions (two-way and three-way)
results1amide1 = factorialANOVA("amide1", "Protein", "Y")

# (2) Testing two-way interactions only, since three-way is insignificant
results2amide1 = factorialANOVA("amide1", "Protein", None, timeXppt, timeXveg,
                                vegXppt)

# (3) Testing for main effects only, since all interactions are insignificant
results3amide1 = factorialANOVA("amide1", "Protein")
# %%
# Purpose: Running factorial ANOVAs on amide2

# (1) Testing for all interactions (two-way and three-way)
results1amide2 = factorialANOVA("amide2", "Protein", "Y")

# (2) Testing for two-way interactions only
results2amide2 = factorialANOVA("amide2", "Protein", None, timeXppt, timeXveg,
                                vegXppt)

# (3) Testing for main effects only
results3amide2 = factorialANOVA("amide2", "Protein")
"""Only precipitation is significant"""
# %%
# Purpose: Reading in the CAZyme domains data and wrangling it

# (1) Reading in the CAZyme domains data
CAZfolder = cwd/"CAZyme metagenomic data"
CAZpath = CAZfolder/"Wrangled CAZyme gene domains.xlsx"
CAZ = pd.read_excel(CAZpath)

# (2) Wrangling the CAZyme domains data by renaming a column and attaching the
# dataframe to the parameters dataframe
CAZ.rename(columns={"Substrate": "Enzyme"}, inplace=True)
parameters = pd.concat([parameters, CAZ])

# %%
# Purpose: Running factorial ANOVAs on hemicellulose CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1hemi = factorialANOVA("Hemicellulose", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2hemi = factorialANOVA("Hemicellulose", "Relative abundance", None,
                              timeXppt, timeXveg, vegXppt)

# (3) Testing for all main effects and only the two-way interaction
# time x vegetation
results3hemi = factorialANOVA("Hemicellulose", "Relative abundance", None,
                              timeXveg)
# %%
# Purpose: Running factorial ANOVAS on lignin CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1lignin = factorialANOVA("Lignin", "Relative abundance", "Y")
"""And there's a significant three-way interaction right off the bat"""
# %%
# Purpose: Running factorial ANOVAs on polysaccharide CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1poly = factorialANOVA("Polysaccharide", "Relative abundance", "Y")
"""And, already, a significant three-way interaction off the bat"""
# %%
# Purpose: Running factorial ANOVAs on oligosaccharide CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1oligo = factorialANOVA("Oligosaccharides", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2oligo = factorialANOVA("Oligosaccharides", "Relative abundance", None,
                               timeXppt, timeXveg, vegXppt)

# (3) Testing for main effects and only the two-way interaction time x veg
results3oligo = factorialANOVA("Oligosaccharides", "Relative abundance", None,
                               timeXveg)
# %%
# Purpose: Running factorial ANOVAs on cell wall CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1cellWall = factorialANOVA("Cell_wall", "Relative abundance", "Y")
"""Another three-way interaction off the bat. Great."""
# %%
# Purpose: Running factorial ANOVAs on Inulin CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1inulin = factorialANOVA("Inulin", "Relative abundance", "Y")
"""Seriously? Another three-way interaction? How am I supposed to interpret
three-way interactions? They're too complicated"""
# %%
# Purpose: Running factorial ANOVAs on starch CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1starch = factorialANOVA("Starch", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2starch = factorialANOVA("Starch", "Relative abundance", None, timeXppt,
                                timeXveg, vegXppt)

# (3) Testing for main effects and only the vegetation x precipitaion two-way
# interaction
results3starch = factorialANOVA("Starch", "Relative abundance", None, vegXppt)
# %%
# Purpose: Running factorial ANOVAs on trehalose CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1trehalose = factorialANOVA("Trehalose", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2trehalose = factorialANOVA("Trehalose", "Relative abundance", None,
                                   timeXppt, timeXveg, vegXppt)

# (3) Testing for main effects only
results3trehalose = factorialANOVA("Trehalose", "Relative abundance")
# %%
# Purpose: Running factorial ANOVAs on cellulose CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1cellulose = factorialANOVA("Cellulose", "Relative abundance", "Y")
"""Another three-way interaction off the bat. Great."""
# %%
# Purpose: Running factorial ANOVAs on pectin CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1pectin = factorialANOVA("Pectin", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2pectin = factorialANOVA("Pectin", "Relative abundance", None, timeXppt,
                                timeXveg, vegXppt)

# (3) Testing for main effects and the two-way interaction of veg x ppt
results3pectin = factorialANOVA("Pectin", "Relative abundance", None, vegXppt)
# %%
# Purpose: Running factorial ANOVAs on glycogen CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1glycogen = factorialANOVA("Glycogen", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2glycogen = factorialANOVA("Glycogen", "Relative abundance", None,
                                  timeXppt, timeXveg, vegXppt)

# (3) Testing for main effects and the time x vegetation interaction only
results3glycogen = factorialANOVA("Glycogen", "Relative abundance", None,
                                  timeXveg)
# %%
# Purpose: Running factorial ANOVAs on peptidoglycan CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1peptido = factorialANOVA("Peptidoglycan", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2peptido = factorialANOVA("Peptidoglycan", "Relative abundance", None,
                                 timeXppt, timeXveg, vegXppt)

# (3) Testing for main effects only
results3peptido = factorialANOVA("Peptidoglycan", "Relative abundance")
# %%
# Purpose: Running factorial ANOVAs on chitin CAZyme domains

# (1) Testing for all interactions (two-way and three-way)
results1chitin = factorialANOVA("Chitin", "Relative abundance", "Y")

# (2) Testing for two-way interactions only
results2chitin = factorialANOVA("Chitin", "Relative abundance", None, timeXppt,
                                timeXveg, vegXppt)

# (3) Testing for all main effects and only the time x veg two-way interaction
results3chitin = factorialANOVA("Chitin", "Relative abundance", None, timeXveg)

# (4) Testing for main effects only
results4chitin = factorialANOVA("Chitin", "Relative abundance")
# %%
# Purpose: Running factorial ANOVAs on total CAZyme domain counts

# (1) Testing for all interactions (two-way and three-way)
results1total = factorialANOVA("Total", "Total CAZyme domains", "Y")
"""And a three-way interaction off the bat. Great."""
# %%
print(datetime.now() - start)
