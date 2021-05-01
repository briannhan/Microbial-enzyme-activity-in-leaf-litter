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

# Reading in parameters file
cwd = Path(os.getcwd())
cwdFilesDirs = os.listdir(cwd)
eeaFolder = cwd/'Enzyme activity data'
eeaFilesDirs = os.listdir(eeaFolder)
paramsPath = eeaFolder/'Parameters.xlsx'
parameters = pd.read_excel(paramsPath)
# %%
# Purpose: Performing the Shapiro-Wilk test

# (1) Creating a dataframe containing only 1 column (enzyme names) that will
# be used as a format for dataframes that hold results
enzymeNames = ["AG", "AP", "BG", "BX", "CBH", "LAP", "NAG", "PPO"]
resultsTemplate = pd.DataFrame({"Enzyme": enzymeNames})


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
    tpParams = paramsDF[paramsDF["Time point"] == timePoint]
    paramsGrouped = tpParams.groupby(["Precip", "Vegetation", "Enzyme",
                                      "Parameter"])
    VmaxResultsDF = resultsTemplate.copy()
    KmResultsDF = resultsTemplate.copy()
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
        if parameter == "Vmax":
            VmaxResultsDF.loc[enzymeInd, column] = resultStr
        elif parameter == "Km":
            KmResultsDF.loc[enzymeInd, column] = resultStr
    return VmaxResultsDF, KmResultsDF


# (3) Perform Shapiro-Wilk test for each enzyme parameter for each time point,
# producing 2 dataframes per time point (1 for each parameter type), resulting
# in 8 dataframes in total. These are test results for when before any data
# transformations had taken place (which I will by "Virgin").
T0VmaxVirgin, T0KmVirgin = shapiroWilk(0, parameters)
T3VmaxVirgin, T3KmVirgin = shapiroWilk(3, parameters)
T5VmaxVirgin, T5KmVirgin = shapiroWilk(5, parameters)
T6VmaxVirgin, T6KmVirgin = shapiroWilk(6, parameters)

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
exportShapiroWilk("No transformations.xlsx", virginDFs)
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
exportShapiroWilk("Natural log transformed.xlsx", lnTransResults)
lnTransParamsPath = eeaFolder/"Parameters - natural log transformed.xlsx"
lnTransParams.to_excel(lnTransParamsPath, index=False)
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
exportShapiroWilk("Log 10 transformed.xlsx", log10TransResults)
log10TransParamsPath = eeaFolder/"Parameters - log 10 transformed.xlsx"
log10TransParams.to_excel(log10TransParamsPath, index=False)
'''Log transforming by base 10 has the same effects as the natural log, with
the exact same results (by asterisks annotation) as indicated by the results
files. So, this data transformation method does improve normality for the most
part, except for T6 parameters as discussed above.'''
