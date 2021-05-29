# -*- coding: utf-8 -*-
"""
Created on Tue May 11 15:32:03 2021

@author: Brian Chung
This script performs Tukey posthoc tests following factorial MANOVA and ANOVA
results. This posthoc analysis, along with the factorial ANOVAs conducted prior
to this analysis, draws heavily from the following links:
https://youtu.be/d_Azlncd-kU
https://www.pythonfordatascience.org/factorial-anova-python/

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
VmaxANOVA = VmaxANOVA[VmaxANOVA.Enzyme != "MANOVA"]
KmANOVA = pd.read_excel(ANOVAresults, "Km").dropna(axis=1)
KmANOVA = KmANOVA[KmANOVA.Enzyme != "MANOVA"]
ANOVAcols = VmaxANOVA.columns.tolist()
mainEffects = ANOVAcols[1:4]
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
# TukeyFile = open(TukeyFilePath, "w")
# for i in range(len(TukeyResults)):
#     row = TukeyResults[i]
#     for j in range(len(row)):
#         value = str(row[j])
#         TukeyFile.write(value)
#         TukeyFile.write(' ')
#     TukeyFile.write('\n')
# TukeyFile.close()

# (4) Writing the Tukey text file to an Excel file
APTukeyResultsName = "AP Tukey.xlsx"
APTukeyPath = statsFolder/"Tukey posthoc"/"AP"/APTukeyResultsName
# APVmaxTukeyDF = pd.read_csv(TukeyFilePath, sep=" ", header=0).dropna(axis=1)
# with ExcelWriter(APTukeyPath) as writer:
#     APVmaxTukeyDF.to_excel(writer, "Vmax", index=False)
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

    ezANOVAind = ANOVAdf[ANOVAdf.Enzyme == ez].index.tolist()[0]
    conditions = (parameters.Enzyme == ez) & (parameters.Parameter == pm)
    paramsOI = parameters[conditions]  # parameters of interest
    interToTest = []
    for interaction in interactions:
        ANOVAvalue = ANOVAdf.loc[ezANOVAind, interaction]
        if ANOVAvalue != "o":
            interToTest.append(interaction)

    sigMainEffect = []
    for mainEffect in mainEffects:
        ANOVAvalue = ANOVAdf.loc[ezANOVAind, mainEffect]
        if ANOVAvalue != "o":
            sigMainEffect.append(mainEffect)

    if "Precipitation" in sigMainEffect:
        precipInd = sigMainEffect.index("Precipitation")
        sigMainEffect[precipInd] = "Precip"

    mainEnotInInter = []
    effectToTest = []
    for mainEffect in sigMainEffect:
        for interaction in interToTest:
            if mainEffect not in interaction:
                mainEnotInInter.append(mainEffect)

    if "Three-way" in interToTest or len(interToTest) >= 2:
        sigMainEffect = []
        mainEnotInInter = []
        effectToTest = []

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
            tukeyTxtName = "{0}, {1}.txt".format(pm, interaction)
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
            tukeyExcelName = "{0}, {1}.xlsx".format(pm, interaction)
            tukeyExcelPath = resultsFolder/tukeyExcelName
            tukeyDF = pd.read_csv(tukeyTxtPath, sep=",", header=0).dropna(1)
            tukeyDF.to_excel(tukeyExcelPath, pm, index=False)
    if len(mainEnotInInter) > 0:
        effectToTest = mainEnotInInter
    elif len(interToTest) == 0 and len(sigMainEffect) > 0:
        effectToTest = sigMainEffect

    if len(effectToTest) > 0:
        # Performing pairwise contrasts Tukey tests for parameters that either
        # have no significant interactions or have significant interactions
        # and 1 main effect that's not part of the interaction
        for mainE in effectToTest:
            groups = paramsOI[mainE]
            tukeyResults = MC(paramsOI.value, groups).tukeyhsd().summary()

            # Writing Tukey results out to a .txt file
            resultsFolder = statsFolder/"Tukey posthoc"/ez
            tukeyTxtName = "{0}, {1}.txt".format(pm, mainE)
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
            tukeyExcelName = "{0}, {1}.xlsx".format(pm, mainE)
            tukeyExcelPath = resultsFolder/tukeyExcelName
            tukeyDF = pd.read_csv(tukeyTxtPath, sep=",", header=0).dropna(1)
            if not os.path.exists(tukeyExcelPath):
                tukeyDF.to_excel(tukeyExcelPath, pm, index=False)
    return


# (2) Performing Tukey tests for all enzymes
# Tukey("AP", "Vmax")
# Tukey("AP", "Km")

# Tukey("BG", "Vmax")
'''Interestingly, BG Vmax shows no significant pairwise differences for time
point despite ANOVA showing time point as a significant main effect. Anyway,
the timePoint excel and text files are deleted manually because of this.
'''

# Tukey("BX", "Km")

# Tukey("CBH", "Vmax")
'''Likewise, CBH Vmax shows no significant pairwise differences for
timePoint x Precipitation despite ANOVA results indicating otherwise, so these
files are deleted manually.

So, I can also look at Precipitation as a simple main effect on CBH Vmax, then.
'''
# Tukey("CBH", "Km")

# Tukey("LAP", "Km")

# Tukey("NAG", "Vmax")
'''And also, NAG Vmax shows no signifcant pairwise differences for
precipitation as a main effect despite ANOVA results indicating otherwise. So,
these are deleted accordingly.
'''
# Tukey("NAG", "Km")

# Tukey("PPO", "Vmax")
'''And PPO Vmax shows no significant pairwise differences for
timePoint x Vegetation as a main effect, so these files are deleted accordingly
'''
# Tukey("PPO", "Km")

'''The files that are deleted manually above show significance at the 0.05
level in ANOVA, but this alpha level is likely too high, resulting in a lack of
significance according to Tukey posthoc.
'''

# (3) Comparing mean BX Km across the 2 vegetation types just to see which one
# has larger Km
BXKm = parameters[(parameters.Enzyme == "BX") & (parameters.Parameter == "Km")]
meanBXKm = BXKm.groupby("Vegetation")["value"].mean()
'''So the way that the tukey results dataframe is set up is that it subtracts
the mean of group 1 from group 2 so that
meandiff = group2 - group1
'''
# %%
# Purpose: Simplifying Tukey test results


# (1) Writing text files of groups in the output Tukey results Excel files to
# facilitate manual annotation of groups with significant differences
def groups(enzyme):
    """
    Generates Excel files of Tukey posthoc groups of all significant
    interactions and significant main effects not part of these interactions
    for both parameters of an enzyme.

    Parameters
    ----------
    enzyme : str
        The enzyme of interest. This specifies the folder that contain Tukey
        posthoc results of a particular enzyme.

    Returns
    -------
    None.

    """
    resultsFolder = statsFolder/"Tukey posthoc"/enzyme
    allResults = os.listdir(resultsFolder)
    if enzyme == "AP":
        '''Some files in the AP Tukey results folder begins with the enzyme
        name while in the other Tukey results folders, no files begin with
        enzyme names. So I'm removing the files that begin with enzyme names
        to facilitate the flow of this function.
        '''
        correctedResults = []
        for file in allResults:
            if "AP" not in file:
                correctedResults.append(file)
        allResults = correctedResults

    # Isolating only Excel files that contain Tukey results to work with
    excelResults = []
    for file in allResults:
        if file.endswith(".xlsx") and "groups" not in file:
            excelResults.append(file)

    # Creating Excel files of significant groups
    for file in excelResults:
        filePath = resultsFolder/file
        results = pd.read_excel(filePath)
        group1 = results.group1.tolist()
        group1 = list(set(group1))  # Dropping duplicates
        group2 = results.group2.tolist()
        group2 = list(set(group2))
        groups = group1 + group2  # Dropping duplicates
        sigGroupsDF = pd.DataFrame({"groups": groups})
        sigGroupsDF = sigGroupsDF.drop_duplicates()  # Dropping duplicates
        rawFileName = file.strip(".xlsx")
        outputName = "{0}, groups.xlsx".format(rawFileName)
        outputPath = resultsFolder/outputName
        if not os.path.exists(outputPath):
            sigGroupsDF.to_excel(outputPath, index=False)
    return


# (2) Creating Excel files of significant groups using the function above
# groups("AP")
# groups("BG")
# groups("BX")
# groups("CBH")
# groups("LAP")
# groups("NAG")
# groups("PPO")
'''And now, some manual annotation of which groups are significant or not.
'''
# %%
# Purpose: Testing Precipitation as a main effect on CBH Vmax

# (1) Reading in new ANOVA results file that's updated with non-significant
# results from Tukey posthoc tests
ANOVApath = statsFolder/'ANOVA, updated with Tukey.xlsx'
ANOVAresults = ExcelFile(ANOVApath)
VmaxANOVA = pd.read_excel(ANOVAresults, "Vmax").dropna(axis=1)
VmaxANOVA = VmaxANOVA[VmaxANOVA.Enzyme != "MANOVA"]
'''VmaxANOVA should easily feed into Tukey()
'''

# (2) Get new Tukey results and Tukey groups
Tukey("CBH", "Vmax")
groups("CBH")
'''For some reason timePoint x Vegetation was significant for ANOVA and is also
significant for Tukey but I didn't have the results file for it, so I didn't
have the boxplot for it as well. But running this chunk just now gives me
results for timePoint x Vegetation and Precipitation as a simple main effect
for CBH Vmax. Interesting.
'''
# %%
print(dt.now() - start)
