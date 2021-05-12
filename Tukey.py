# -*- coding: utf-8 -*-
"""
Created on Tue May 11 15:32:03 2021

@author: Brian Chung
This script performs Tukey posthoc tests following factorial MANOVA and ANOVA
results. This posthoc analysis, along with the factorial ANOVAs conducted prior
to this analysis, draws heavily from the following YouTube video:
https://youtu.be/d_Azlncd-kU

For significant interactions, simple main effects will be investigated. If
no interactions are significant, then pairwise contrasts between the main
effects will be investigated. These differences will all be investigated using
the Tukey Honest Significant Difference (HSD) posthoc test.
"""

# (1) Importing necessary libraries to load in and wrangle data, time how long
# it takes to run this script (for my curiosity really), and to conduct
# Tukey posthoc tests.
import pandas as pd
from pandas import ExcelFile
from pandas import ExcelWriter
import os
from datetime import datetime as dt
from statsmodels.sandbox.stats.multicomp import MultiComparison as MC
from pathlib import Path

# (2) Reading in ANOVA results
start = dt.now()
cwd = Path(os.getcwd())
cwdDirsFiles = os.listdir(cwd)
statsFolder = cwd/'Statistical analyses'
statsContents = os.listdir(statsFolder)
ANOVApath = statsFolder/'ANOVA results.xlsx'
ANOVAresults = ExcelFile(ANOVApath)
VmaxANOVA = pd.read_excel(ANOVAresults, "Vmax").dropna(axis=1)
KmANOVA = pd.read_excel(ANOVAresults, "Km").dropna(axis=1)
ANOVAcols = VmaxANOVA.columns.tolist()
mainEffects = ANOVAcols[1:4]
for i in range(len(mainEffects)):
    mainEffect = mainEffects[i]
    if mainEffect == "Precipitation":
        # Changes this name to have the same name as a column in the parameters
        # dataframe below in order to more easily subset the parameters
        # dataframe
        mainEffects[i] = "Precip"
interactions = ANOVAcols[4:]
twoWay = interactions[:-1]

# (3) Reading in log 10 transformed parameters
activityFolder = cwd/'Enzyme activity data'
activityContents = os.listdir(activityFolder)
paramsPath = activityFolder/'Parameters - log 10 transformed.xlsx'
parameters = pd.read_excel(paramsPath).drop("Transformation", 1)
parameters.timePoint = parameters.timePoint.astype(str)
# %%
# Purpose: Tukey posthoc on AP Vmax as a proof of concept before I can
# generalize Tukey posthoc into a function or series of functions

# (1) Isolating AP Vmax from main parameters dataframe
AP_Vmax_ANOVA = VmaxANOVA[VmaxANOVA.Enzyme == "AP"]
APconditions = (parameters.Enzyme == "AP") & (parameters.Parameter == "Vmax")
VmaxAP = parameters[APconditions]

# (2) Now I'm gonna test for timePoint as a main effect
timePoints = VmaxAP.timePoint
comparisons = MC(VmaxAP.value, timePoints)
TukeyResults = comparisons.tukeyhsd().summary()

# (3) Writing Tukey results to a text file
TukeyFilePath = statsFolder/"Tukey posthoc"/"AP"/"AP Vmax Tukey.txt"
TukeyFile = open(TukeyFilePath, "w")
for i in range(len(TukeyResults)):
    row = TukeyResults[i]
    for j in range(len(row)):
        value = str(row[j])
        TukeyFile.write(value)
        TukeyFile.write(' ')
    TukeyFile.write('\n')
TukeyFile.close()

# (4) Writing the Tukey text file to an Excel file
APTukeyResultsName = "AP Tukey.xlsx"
APTukeyPath = statsFolder/"Tukey posthoc"/"AP"/APTukeyResultsName
APVmaxTukeyDF = pd.read_csv(TukeyFilePath, sep=" ", header=0).dropna(axis=1)
with ExcelWriter(APTukeyPath) as writer:
    APVmaxTukeyDF.to_excel(writer, "Vmax", index=False)
# %%
# Purpose: Performing Tukey tests for all enzymes (except for AG, which has
# no significant interactions or main effects for both Vmax and Km)

# (1) Abstracting Tukey posthoc testing into a single function or a few
# functions


def Tukey(ez, pm):
    """
    Performs Tukey posthoc tests on a parameter of a particular enzyme and
    reads the Tukey results out to a text file and an Excel file.

    Parameters
    ----------
    ez : str
        The enzyme of interest, used to isolate parameter values specific to
        this enzyme as well as ANOVA results specific to this enzyme. Also used
        to name results files.
    pm : str
        Parameter of interest (either Vmax or Km). Used to indicate the name of
        the sheet in the Excel file that will contain the Tukey results.

    Returns
    -------
    None.

    """
    if pm == "Vmax":
        ANOVAdf = VmaxANOVA
    elif pm == "Km":
        ANOVAdf = KmANOVA

    enzymeANOVA = ANOVAdf[ANOVAdf.Enzyme == ez]
    conditions = (parameters.Enzyme == ez) & (parameters.Parameter == pm)
    paramsOI = parameters[conditions]  # parameters of interest
    interToTest = []
    for interaction in interactions:
        ANOVAvalue = str(enzymeANOVA[interaction])
        if ANOVAvalue != "o":
            interToTest.append(interaction)

    mainEffectToTest = []
    if len(interToTest) == 0:
        for mainEffect in mainEffects:
            ANOVAvalue = str(enzymeANOVA[mainEffect])
            if ANOVAvalue != "o":
                mainEffectToTest.append(mainEffect)

    if len(interToTest) > 0:  # Perform simple main effects Tukey test
        # print(interToTest)
        for interaction in interToTest:
            if interaction == "Three-way":
                groups = (paramsOI.timePoint + " x " + paramsOI.Precip + " x "
                          + paramsOI.Vegetation)
            elif interaction in twoWay:
                mainEsplit = interaction.split(' x ')
                for i in range(2):
                    if mainEsplit[i] == "Precipitation":
                        mainEsplit[i] = "Precip"
                groups = paramsOI[mainEsplit[0]]+" x "+paramsOI[mainEsplit[1]]
            tukeyResults = MC(paramsOI.value, groups).tukeyhsd().summary()

            # Writing Tukey results out to a .txt file
            resultsFolder = statsFolder/"Tukey posthoc"/ez
            tukeyTxtName = interaction + ".txt"
            tukeyTxtPath = resultsFolder/tukeyTxtName
            txtFile = open(tukeyTxtPath, "w")
            resultsLen = len(tukeyResults)
            for i in range(resultsLen):
                row = tukeyResults[i]
                for value in row:
                    value = str(value)
                    txtFile.write(value)
                    txtFile.write(",")
                txtFile.write("\n")
            txtFile.close()

            # Writing Tukey results out to an Excel file
            tukeyExcelName = interaction + ".xlsx"
            tukeyExcelPath = resultsFolder/tukeyExcelName
            tukeyDF = pd.read_csv(tukeyTxtPath, sep=",", header=0)
            with ExcelWriter(tukeyExcelPath) as writer:
                tukeyDF.to_excel(writer, pm)
    elif len(interToTest) == 0 and len(mainEffectToTest) > 0:
        # print(mainEffectToTest)
        # Performing pairwise contrasts Tukey tests
        for mainE in mainEffects:
            groups = paramsOI[mainE]
            tukeyResults = MC(paramsOI.value, groups).tukeyhsd().summary()

            # Writing Tukey results out to a .txt file
            resultsFolder = statsFolder/"Tukey posthoc"/ez
            tukeyTxtName = mainE + ".txt"
            tukeyTxtPath = resultsFolder/tukeyTxtName
            txtFile = open(tukeyTxtPath, "w")
            resultsLen = len(tukeyResults)
            for i in range(resultsLen):
                row = tukeyResults[i]
                for value in row:
                    value = str(value)
                    txtFile.write(value)
                    txtFile.write(" ")
                txtFile.write("\n")
            txtFile.close()

            # Writing Tukey results out to an Excel file
            tukeyExcelName = interaction + ".xlsx"
            tukeyExcelPath = resultsFolder/tukeyExcelName
            tukeyDF = pd.read_csv(tukeyTxtPath, sep=",", header=0).dropna(1)
            with ExcelWriter(tukeyExcelPath) as writer:
                tukeyDF.to_excel(writer, pm, index=False)
    return


# (2) Performing Tukey tests for all enzymes
# print("AP Vmax interactions/main effects to test:")
Tukey("AP", "Vmax")
# print('\n')

# print("AP Km interactions/main effects to test")
Tukey("AP", "Km")
# print('\n')

# print("BG Vmax interactions/main effects to test")
Tukey("BG", "Vmax")
# print('\n')
# %%
print(dt.now() - start)
