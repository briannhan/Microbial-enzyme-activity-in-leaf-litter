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
def shapiroWilk(timePoint):
    """
    Performs the Shapiro-Wilk test for normality on enzyme parameters from a
    particular time point, outputting 2 dataframes, 1 for each parameter.

    Parameters
    ----------
    timePoint : int
        A number that's unique to a particular set of parameters, with the
        parameters having been estimated from data from this time point.
    VmaxResultsDF : Pandas dataframe
        A dataframe that will hold the Vmax results for this particular time
        point.
    KmResultsDF : Pandas dataframe
        A dataframe that will hold the Km results for this particular time
        point.

    Returns
    -------
    2 dataframes of test results, 1 for each parameter.

    """
    tpParams = parameters[parameters["Time point"] == timePoint]
    paramsGrouped = tpParams.groupby(["Precip", "Vegetation", "Enzyme"])
    VmaxResultsDF = resultsTemplate.copy()
    KmResultsDF = resultsTemplate.copy()
    for name, group in paramsGrouped:
        column = "{0}, {1}".format(name[0], name[1])
        Vmax = group["Vmax"]
        testVmax = stats.shapiro(Vmax)
        VmaxP = testVmax.pvalue
        if VmaxP < 0.05 and VmaxP >= 0.01:
            VmaxResult = "*"
        elif VmaxP < 0.01 and VmaxP >= 0.001:
            VmaxResult = "**"
        elif VmaxP < 0.001:
            VmaxResult = "***"
        elif VmaxP >= 0.05:
            VmaxResult = "o"
        enzymeInd = VmaxResultsDF[VmaxResultsDF["Enzyme"] == name[2]].index
        VmaxResultsDF.loc[enzymeInd, column] = VmaxResult

        Km = group["Km"]
        testKm = stats.shapiro(Km)
        KmP = testKm.pvalue
        if KmP < 0.05 and KmP >= 0.01:
            KmResult = "*"
        elif KmP < 0.01 and KmP >= 0.001:
            KmResult = "**"
        elif KmP < 0.001:
            KmResult = "***"
        elif KmP >= 0.05:
            KmResult = "o"
        KmResultsDF.loc[enzymeInd, column] = KmResult
    return VmaxResultsDF, KmResultsDF


# (3) Perform Shapiro-Wilk test for each enzyme parameter for each time point,
# producing 2 dataframes per time point (1 for each parameter type), resulting
# in 8 dataframes in total. These are test results for when before any data
# transformations had taken place (which I will by "Virgin").
T0VmaxVirgin, T0KmVirgin = shapiroWilk(0)
T3VmaxVirgin, T3KmVirgin = shapiroWilk(3)
T5VmaxVirgin, T5KmVirgin = shapiroWilk(5)
T6VmaxVirgin, T6KmVirgin = shapiroWilk(6)

# (4) Logistics for exporting the virgin results, including declaring the path
# to save the results and writing out the names for each sheet of the results
statsResultsFolder = cwd/'Statistical analyses'
statsResultsList = os.listdir(statsResultsFolder)
normTestFolder = statsResultsFolder/'Normality testing'
virginResults = normTestFolder/"No transformations.xlsx"
virginDFs = [T0VmaxVirgin, T0KmVirgin, T3VmaxVirgin, T3KmVirgin, T5VmaxVirgin,
             T5KmVirgin, T6VmaxVirgin, T6KmVirgin]
resultsNames = ["T0 Vmax", "T0 Km", "T3 Vmax", "T3 Km", "T5 Vmax", "T5 Km",
                "T6 Vmax", "T6 Km"]


# (5) Writing a function to export the Shapiro-Wilk test results
def exportShapiroWilk(fileName, resultsDFs, resultsNames):
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
            
