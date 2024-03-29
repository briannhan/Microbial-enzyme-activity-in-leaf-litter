# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 14:48:06 2021

@author: Brian Chung
This script performs the Shapiro-Wilk test for normality on the parameters from
each treatment combination from each timepoint. For example, it will test for
the normality of Vmax from grassland x drought in T0, and then CSS x drought
for T0, and so on, do forth for T0, and then T3, T5, and T6. It will produce
Excel files of Shapiro-Wilk test results, and then I will use these results to
transform the parameters until the parameters are not significantly different
according to the Shapiro-Wilk test.
"""
import pandas as pd
from scipy import stats
import os
from pathlib import Path
from pandas import ExcelWriter
import numpy as np
from matplotlib import pyplot as py

# Reading in parameters file
cwd = Path(os.getcwd())
cwdFilesDirs = os.listdir(cwd)
eeaFolder = cwd/'Enzyme activity data'
eeaFilesDirs = os.listdir(eeaFolder)
paramsPath = eeaFolder/'Parameters.xlsx'
parameters = pd.read_excel(paramsPath)

# Reading in wrangled litter chemistry data
litterChemFolder = cwd/"Litter chemistry"
litterChemPath = litterChemFolder/"Carbohydrates and Proteins FTIR.xlsx"
litterChem = pd.read_excel(litterChemPath)

# Renaming columns in the litter chemistry data to apply the existing
# Shapiro-Wilk code to it
litterChem["Parameter"] = "FTIR spectral area"
litterChem.rename(columns={"functionalGroup": "Enzyme"}, inplace=True)
litterChem["timePoint"] = litterChem["timePoint"].astype(int)

# Reading in the CAZymes gene domain relative abundance data
CAZfolder = cwd/"CAZyme metagenomic data"
CAZpath = CAZfolder/"Wrangled CAZyme gene domains.xlsx"
CAZ = pd.read_excel(CAZpath)

# Renaming columns in the CAZymes gene domain data
CAZ.rename(columns={"Substrate": "Enzyme"}, inplace=True)
CAZ["timePoint"] = CAZ["timePoint"].astype(int)
# %%
# Purpose: Performing the Shapiro-Wilk test

# (1) Creating a dataframe containing only 1 column (enzyme names, litter
# chemistry functional groups, or CAZyme substrates) that will
# be used as a format for dataframes that hold results
enzymeNames = ["AG", "AP", "BG", "BX", "CBH", "LAP", "NAG", "PPO"]
resultsTemplate = pd.DataFrame({"Enzyme": enzymeNames})
functionalGroups = list(set(litterChem.Enzyme.tolist()))
resultsTemplateLC = pd.DataFrame({"functionalGroup": functionalGroups})
substrates = list(set(CAZ.Enzyme.tolist()))
resultsTemplateCAZ = pd.DataFrame({"Substrate": substrates})
CAZparameters = set(CAZ.Parameter.tolist())
# LC = litter chemistry


# (2) Creating a function to perform the Shapiro-Wilk test and output 2
# dataframes of test results -- 1 for each particular enzyme parameter --
# from 1 time point
def shapiroWilk(timePoint, paramsDF):
    """
    Performs the Shapiro-Wilk test for normality on enzyme parameters from a
    particular time point, outputting 2 dataframes, 1 for each parameter.

    Parameters
    ----------
    timePoint : int
        A number that's unique to a particular set of parameters, with the
        parameters having been estimated from data from this time point.
    paramsDF : Pandas dataframe
        Dataframe of parameter values for all time points.

    Returns
    -------
    2 dataframes of test results, 1 for each parameter.

    """
    tpParams = paramsDF[paramsDF["timePoint"] == timePoint]
    paramsGrouped = tpParams.groupby(["Precip", "Vegetation", "Enzyme",
                                      "Parameter"])
    VmaxResultsDF = resultsTemplate.copy()
    KmResultsDF = resultsTemplate.copy()
    lcResultsDF = resultsTemplateLC.copy()  # lc = litter chemistry
    CAZresultsDF = resultsTemplateCAZ.copy()
    parametersToTest = set(tpParams.Parameter.tolist())
    for name, group in paramsGrouped:
        column = "{0}, {1}".format(name[0], name[1])
        enzyme = name[2]
        parameter = name[-1]
        values = group["value"]
        test = stats.shapiro(values)
        p = test.pvalue
        if p < 0.05 and p >= 0.01:
            resultStr = "*"
        elif p < 0.01 and p >= 0.001:
            resultStr = "**"
        elif p < 0.001:
            resultStr = "***"
        elif p >= 0.05:
            resultStr = "o"
        enzymeInd = VmaxResultsDF[VmaxResultsDF["Enzyme"] == enzyme].index
        fgInd = lcResultsDF[lcResultsDF.functionalGroup == enzyme].index
        CAZind = CAZresultsDF[CAZresultsDF.Substrate == enzyme].index
        if parameter == "Vmax":
            VmaxResultsDF.loc[enzymeInd, column] = resultStr
        elif parameter == "Km":
            KmResultsDF.loc[enzymeInd, column] = resultStr
        elif parameter == "FTIR spectral area":
            lcResultsDF.loc[fgInd, column] = resultStr
        elif parameter in CAZparameters:
            CAZresultsDF.loc[CAZind, column] = resultStr
    if "Vmax" in parametersToTest or "Km" in parametersToTest:
        return VmaxResultsDF, KmResultsDF
    elif "FTIR spectral area" in parametersToTest:
        return lcResultsDF
    else:
        return CAZresultsDF


# (3) Perform Shapiro-Wilk test for each enzyme parameter for each time point,
# producing 2 dataframes per time point (1 for each parameter type), resulting
# in 8 dataframes in total. These are test results for when before any data
# transformations had taken place (which I will by "Virgin").
T0VmaxVirgin, T0KmVirgin = shapiroWilk(0, parameters)
T3VmaxVirgin, T3KmVirgin = shapiroWilk(3, parameters)
T5VmaxVirgin, T5KmVirgin = shapiroWilk(5, parameters)
T6VmaxVirgin, T6KmVirgin = shapiroWilk(6, parameters)
# %%
# (4) Logistics for exporting the virgin results, including declaring the path
# to save the results and writing out the names for each sheet of the results
statsResultsFolder = cwd/'Statistical analyses'
statsResultsList = os.listdir(statsResultsFolder)
normTestFolder = statsResultsFolder/'Normality testing'
virginDFs = [T0VmaxVirgin, T0KmVirgin, T3VmaxVirgin, T3KmVirgin, T5VmaxVirgin,
             T5KmVirgin, T6VmaxVirgin, T6KmVirgin]
resultsNames = ["T0 Vmax", "T0 Km", "T3 Vmax", "T3 Km", "T5 Vmax", "T5 Km",
                "T6 Vmax", "T6 Km"]


# (5) Writing a function to export the Shapiro-Wilk test results
def exportShapiroWilk(fileName, resultsDFs):
    """
    Exports an Excel file of Shapiro-Wilk test results. This function is
    intended to be used prior to any data transformation and after each round
    of data transformations.

    Parameters
    ----------
    fileName : str
        The name of the Excel file. Must have the ".xlsx" extension at the
        end.
    resultsDFs : list of Pandas dataframes
        A list of dataframes of test results either before any data
        transformation or after each round of data transformation. Will be
        exported as a new Excel file, with each dataframe having a sheet to its
        own.

    Returns
    -------
    None.

    """
    filePath = normTestFolder/fileName
    with ExcelWriter(filePath) as writer:
        for i in range(8):
            results = resultsDFs[i]
            sheet = resultsNames[i]
            results.to_excel(writer, sheet, index=False)
    return


# (6) Export virgin Shapiro-Wilk results
# exportShapiroWilk("No transformations.xlsx", virginDFs)
# %%
# Purpose: Transforming parameters by taking the natural log of each parameter
'''Well, I'm thinking that, if I transform the data, I would need to transform
all of the data, not just each treatment combination. Because if for a
particular enzyme parameter, it has 1 non-normal treatment combination, and I
transform only that combination instead of the other 3 combinations as well,
then this transformation of only 1 combination can result in false positives
in my ANOVAs due to the transformation causing the non-normal treatment
combination to be different from the other treatment combinations rather than
all treatment combinations, including the non-normal one, being different to
begin with.

Also, how would I perform the t-test between 2 treatment combinations if 1
of them was transformed and the other wasn't? I can't. So, either transform
all or not transform at all.'''

# (1) Log-transforming all parameters
lnTransParams = parameters.copy()
lnTransParams["value"] = np.log(lnTransParams["value"])
lnTransParams["Transformation"] = "Natural log"

# (2) Performing Shapiro-Wilk test on log-transformed parameters
T0VmaxLn, T0KmLn = shapiroWilk(0, lnTransParams)
T3VmaxLn, T3KmLn = shapiroWilk(3, lnTransParams)
T5VmaxLn, T5KmLn = shapiroWilk(5, lnTransParams)
T6VmaxLn, T6KmLn = shapiroWilk(6, lnTransParams)
'''Well, guess that log-transforming the parameters still leave quite a few
as non-normal. Although for the most part, this process does improve normality,
though T6 parameters show mixed results in how normality is affected, with some
treatment combinations seeing improvements in normality while others seeing
normality worsened.'''

# (3) Exporting log-transformed parameters & Shapiro-Wilk test results
lnTransResults = [T0VmaxLn, T0KmLn, T3VmaxLn, T3KmLn, T5VmaxLn, T5KmLn,
                  T6VmaxLn, T6KmLn]
# exportShapiroWilk("Natural log transformed.xlsx", lnTransResults)
lnTransParamsPath = eeaFolder/"Parameters - natural log transformed.xlsx"
# lnTransParams.to_excel(lnTransParamsPath, index=False)
# %%
# Purpose: Transforming parameters by taking the log base 10 of each parameter

# (1) Log-base-10-transforming all parameters
log10TransParams = parameters.copy()
log10TransParams["value"] = np.log10(log10TransParams["value"])
log10TransParams["Transformation"] = "Log10"

# (2) Performing Shapiro-Wilk test on log-transformed parameters
T0VmaxLog10, T0KmLog10 = shapiroWilk(0, log10TransParams)
T3VmaxLog10, T3KmLog10 = shapiroWilk(3, log10TransParams)
T5VmaxLog10, T5KmLog10 = shapiroWilk(5, log10TransParams)
T6VmaxLog10, T6KmLog10 = shapiroWilk(6, log10TransParams)

# (3) Exporting log-transformed parameters & Shapiro-Wilk test results
log10TransResults = [T0VmaxLog10, T0KmLog10, T3VmaxLog10, T3KmLog10,
                     T5VmaxLog10, T5KmLog10, T6VmaxLog10, T6KmLog10]
# exportShapiroWilk("Log 10 transformed.xlsx", log10TransResults)
log10TransParamsPath = eeaFolder/"Parameters - log 10 transformed.xlsx"
# log10TransParams.to_excel(log10TransParamsPath, index=False)
'''Log transforming by base 10 has the same effects as the natural log, with
the exact same results (by asterisks annotation) as indicated by the results
files. So, this data transformation method does improve normality for the most
part, except for T6 parameters as discussed above.'''
# %%
# Purpose: Performing Shapiro-Wilk on untransformed litter chemistry data

# (1) Performing the Shapiro-Wilk test
T0litterChemVirgin = shapiroWilk(0, litterChem)
T3litterChemVirgin = shapiroWilk(3, litterChem)
T5litterChemVirgin = shapiroWilk(5, litterChem)
T6litterChemVirgin = shapiroWilk(6, litterChem)

# (2) Exporting the Shapiro-Wilk test results
noTransformPath = normTestFolder/"Litter chemistry - no transformations.xlsx"
if os.path.exists(noTransformPath) is False:
    with ExcelWriter(noTransformPath) as writer:
        T0litterChemVirgin.to_excel(writer, "T0", index=False)
        T3litterChemVirgin.to_excel(writer, "T3", index=False)
        T5litterChemVirgin.to_excel(writer, "T5", index=False)
        T6litterChemVirgin.to_excel(writer, "T6", index=False)
"""The functional groups that exhibit non-normal behavior are:
lipid: T3 ambient grassland
carboEster1: T0 reduced CSS, T5 reduced CSS
carboEster2: T3 ambient grassland, T6 reduced CSS
total carbohydrate ester: T5 ambient CSS, T5 reduced CSS
amide2: T3 reduced CSS, T5 reduced CSS, T6 ambient CSS
total amide: T5 reduced CSS, T6 reduced CSS
carbohydrate C-O stretching: T6 ambient CSS, T6 reduced grassland

This constitutes almost all functional groups. So I'll apply the Box-Cox
transformation to these functional groups, then re-analyze them with factorial
ANOVAs and Tukey's. Pain peko
"""
# %%
# Purpose: Performing the Box-Cox transformation on leaf litter chemistry
# functional groups that exhibit non-normal behavior

# (1) Performing the Box-Cox transformation on these functional groups
fg2trans = ["lipid", "carboEster1", "carboEster2", "carboEster", "amide2",
            "amide", "C_O_stretching"]
# list of functional groups to transform using Box-Cox
lcBoxCox = litterChem.copy()  # litter chemistry data that will be transformed
for fg in fg2trans:
    untransformed = lcBoxCox.loc[lcBoxCox.Enzyme == fg, "value"]
    transformed, exponent = stats.boxcox(untransformed)
    lcBoxCox.loc[lcBoxCox.Enzyme == fg, "value"] = transformed
    lcBoxCox.loc[lcBoxCox.Enzyme == fg, "exponent"] = exponent

# (2) Performing Shapiro-Wilk on the transformed litter chemistry data
T0litterChemBoxCox = shapiroWilk(0, lcBoxCox)
T3litterChemBoxCox = shapiroWilk(3, lcBoxCox)
T5litterChemBoxCox = shapiroWilk(5, lcBoxCox)
T6litterChemBoxCox = shapiroWilk(6, lcBoxCox)

# (3) Exporting the Shapiro-Wilk test results of the Box-Cox transformed litter
# chemistry data
lcBoxCoxPath = normTestFolder/"Litter chemistry - Box-Cox results.xlsx"
if os.path.exists(lcBoxCoxPath) is False:
    with ExcelWriter(lcBoxCoxPath) as excelFile:
        T0litterChemBoxCox.to_excel(excelFile, "T0", index=False)
        T3litterChemBoxCox.to_excel(excelFile, "T3", index=False)
        T5litterChemBoxCox.to_excel(excelFile, "T5", index=False)
        T6litterChemBoxCox.to_excel(excelFile, "T6", index=False)
"""Ok well, overall, the Box-Cox transformation doesn't seem to improve
normality that much for each individual treatment group. So, I won't be
re-analyzing the litter chemistry data."""
# %%
# Purpose: Testing the normality of untransformed CAZ data

# (1) Performing Shapiro-Wilk on the untransformed CAZyme data
T0CAZvirgin = shapiroWilk(0, CAZ)
T3CAZvirgin = shapiroWilk(3, CAZ)
T5CAZvirgin = shapiroWilk(5, CAZ)
T6CAZvirgin = shapiroWilk(6, CAZ)

# (2) Exporting the Shapiro-Wilk test results
CAZnoTransformPath = normTestFolder/"CAZyme domains - No transformations.xlsx"
if os.path.exists(CAZnoTransformPath) is False:
    with ExcelWriter(CAZnoTransformPath) as CAZresultsFile:
        T0CAZvirgin.to_excel(CAZresultsFile, "T0", index=False)
        T3CAZvirgin.to_excel(CAZresultsFile, "T3", index=False)
        T5CAZvirgin.to_excel(CAZresultsFile, "T5", index=False)
        T6CAZvirgin.to_excel(CAZresultsFile, "T6", index=False)
"""The substrates whose domain relative abundance show deviations from
normality are

Inulin: T3 Ambient CSS, T6 Reduced Grassland
Hemicellulose: T5 Reduced Grassland
Trehalose: T3 Ambient CSS
Peptidoglycan: T6 Ambient CSS
Pectin: T6 Reduced Grassland

Time to perform Box-Cox transformations on these substrates
"""
# %%
# Purpose: Performing Box-Cox transformations on the substrates whose domain
# relative abundance shows deviations from normality and testing them again
# after their respective transformations

# (1) Performing Box-Cox transformations on the above substrates
substrateTrans = ["Inulin", "Hemicellulose", "Trehalose", "Peptidoglycan",
                  "Pectin"]
CAZtransformState = "Not transformed"
"""While performing a Box-Cox transformation, there's a value that's exactly
0: T032X Inulin at time point 0, CSS Ambient. I'm setting that value to
0.00001"""
CAZ.loc[CAZ["value"] == 0, "value"] = 0.00001
if CAZtransformState == "Not transformed":
    for substrate in substrateTrans:
        untransformed = CAZ.loc[CAZ.Enzyme == substrate, "value"]
        transformed, exponent = stats.boxcox(untransformed)
        CAZ.loc[CAZ.Enzyme == substrate, "value"] = transformed
        CAZ.loc[CAZ.Enzyme == substrate, "exponent"] = exponent
    CAZtransformState = "Transformed"

# (2) Performing Shapiro-Wilk on the Box-Cox transformed substrates
T0CAZboxCox = shapiroWilk(0, CAZ)
T3CAZboxCox = shapiroWilk(3, CAZ)
T5CAZboxCox = shapiroWilk(5, CAZ)
T6CAZboxCox = shapiroWilk(6, CAZ)

# (3) Exporting the Shapiro-Wilk test results
CAZboxCoxPath = normTestFolder/"CAZyme domains - Box-Cox.xlsx"
if os.path.exists(CAZboxCoxPath) is False:
    with ExcelWriter(CAZboxCoxPath) as CAZboxCoxResults:
        T0CAZboxCox.to_excel(CAZboxCoxResults, "T0", index=False)
        T3CAZboxCox.to_excel(CAZboxCoxResults, "T3", index=False)
        T5CAZboxCox.to_excel(CAZboxCoxResults, "T5", index=False)
        T6CAZboxCox.to_excel(CAZboxCoxResults, "T6", index=False)
"""Ok well, the transformation didn't improve normality for these substrates.
So, the data is already normal prior to the transformation. Let's do factorial
ANOVAs and Tukey's, now."""

# (4) Testing a Q-Q plot as a way to evaluate normality. I'll make a Q-Q plot
# of the total CAZyme domains.
CAZtotal = CAZ.loc[CAZ.Enzyme == "Total", "value"]
totalQQplot = stats.probplot(CAZtotal, plot=py)
"""Doesn't look quite normal to me, but Shapiro-Wilk states that this is
normal. Nande? Could it be because, I'm plotting a probability plot on ALL of
the 'total' values across all 3 independent variables, so I'm effectively
plotting multiple populations onto a plot for a single population, making it
look like they're not normal because they're from multiple populations?

Well, this explanation makes sense, actually. I'm plotting multiple populations
onto a single plot."""
