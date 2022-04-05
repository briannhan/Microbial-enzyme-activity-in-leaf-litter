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
import cld

# (2) Reading in ANOVA results for enzyme activity (Vmax & Km), litter
# chemistry, and CAZyme domain relative abundance
start = dt.now()
cwd = Path(os.getcwd())
cwdDirsFiles = os.listdir(cwd)
statsFolder = cwd/'Statistical analyses'
statsContents = os.listdir(statsFolder)
ANOVApath = statsFolder/'ANOVA results.xlsx'
ANOVAresults = ExcelFile(ANOVApath)
# Enzyme activity ANOVA results
VmaxANOVA = pd.read_excel(ANOVAresults, "Vmax").dropna(axis=1)
VmaxANOVA = VmaxANOVA[VmaxANOVA.Enzyme != "MANOVA"]
KmANOVA = pd.read_excel(ANOVAresults, "Km").dropna(axis=1)
KmANOVA = KmANOVA[KmANOVA.Enzyme != "MANOVA"]
# ANOVA main effects and interactions
ANOVAcols = VmaxANOVA.columns.tolist()
mainEffects = ANOVAcols[1:4]
interactions = ANOVAcols[4:]
twoWay = interactions[:-1]
# Litter chemistry ANOVA results
litterChemANOVA = pd.read_excel(ANOVAresults, "litterChemistry").dropna(1)
functionalGroups = litterChemANOVA.functionalGroup.tolist()
carboFG = functionalGroups[:3]
miscFG = functionalGroups[3:5]
litterChemANOVA.rename(columns={"functionalGroup": "Enzyme"}, inplace=True)
# CAZyme domains ANOVA results
CAZymeANOVA = pd.read_excel(ANOVAresults, "CAZyme domains").dropna(1)
substrates = CAZymeANOVA.Substrate.tolist()
CAZymeANOVA.rename(columns={"Substrate": "Enzyme"}, inplace=True)

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
    elif pm in functionalGroups:
        ANOVAdf = litterChemANOVA
    elif pm == "Relative abundance" or pm == "Total CAZyme domains":
        ANOVAdf = CAZymeANOVA

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

    if len(interToTest) > 0:  # Perform Tukey test on interactions
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
            if ez != "Total":
                # This if-elif statement is for dealing with total CAZyme
                # domain counts. The name (Total) taken from the
                # ANOVA results and the wrangled data does not match the name
                # of the folder I created (Total CAZyme)
                resultsFolder = statsFolder/"Tukey posthoc"/ez
            elif ez == "Total":
                resultsFolder = statsFolder/"Tukey posthoc"/"Total CAZyme"
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
            if os.path.exists(tukeyExcelPath) is False:
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
        rawFileName = file.rstrip(".xlsx")
        outputName = "{0}, groups.xlsx".format(rawFileName)
        outputPath = resultsFolder/outputName
        if not os.path.exists(outputPath):
            sigGroupsDF.to_excel(outputPath, index=False)
        # cldDF = cld.main(results, 0.05)
        # cldName = "{0}, groups, annotated.xlsx".format(rawFileName)
        # cldPath = resultsFolder/cldName
        # cldDF.to_excel(cldPath, index=False)
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
# ANOVApath = statsFolder/'ANOVA, updated with Tukey.xlsx'
# ANOVAresults = ExcelFile(ANOVApath)
# VmaxANOVA = pd.read_excel(ANOVAresults, "Vmax").dropna(axis=1)
# VmaxANOVA = VmaxANOVA[VmaxANOVA.Enzyme != "MANOVA"]
'''VmaxANOVA should easily feed into Tukey()
'''

# (2) Get new Tukey results and Tukey groups
# Tukey("CBH", "Vmax")
# groups("CBH")
'''For some reason timePoint x Vegetation was significant for ANOVA and is also
significant for Tukey but I didn't have the results file for it, so I didn't
have the boxplot for it as well. But running this chunk just now gives me
results for timePoint x Vegetation and Precipitation as a simple main effect
for CBH Vmax. Interesting.

Turns out, Precipitation is insignificant by Tukey posthoc despite really high
significance as a main effect by ANOVA. So, Tukey posthoc results for
Precipitation are deleted.
'''
# %%
# Purpose: Testing Vegetation as a main effect on PPO Vmax. This is because I
# did not find much of a significant difference when testing for a three-way
# interaction, despite ANOVA stating that there is a significant three-way
# interaction. There were also no other significant interactions or main
# effects aside from Vegetation.

# (1) Performing Tukey's on PPO Vmax and exporting the results
VmaxPPO = parameters[(parameters.Parameter == "Vmax")
                     & (parameters.Enzyme == "PPO")]
vegPPOresults = MC(VmaxPPO.value, VmaxPPO.Vegetation).tukeyhsd().summary()

# PPOfolder = statsFolder/"Tukey posthoc"/"PPO"
# PPOvegTukeyTxtName = "{0}, {1}.txt".format("Vmax", "Vegetation")
# PPOvegTukeyTxtPath = PPOfolder/PPOvegTukeyTxtName
# PPOtxtFile = open(PPOvegTukeyTxtPath, "w")
# PPOresultsRows = len(vegPPOresults)
# for i in range(PPOresultsRows):
#     row = vegPPOresults[i]
#     for value in row:
#         value = str(value)
#         PPOtxtFile.write(value)
#         PPOtxtFile.write(",")
#     PPOtxtFile.write("\n")
# PPOtxtFile.close()

# Writing Tukey results out to an Excel file
# PPOvegTukeyExcelName = "{0}, {1}.xlsx".format("Vmax", "Vegetation")
# PPOvegTukeyExcelPath = PPOfolder/PPOvegTukeyExcelName
# PPOvegTukeyDF = pd.read_csv(PPOvegTukeyTxtPath, sep=",", header=0).dropna(1)
# PPOvegTukeyDF.to_excel(PPOvegTukeyExcelPath, "Vmax", index=False)

# (2) Exporting groups
# groups("PPO")
# %%
# Purpose: Writing a function that abstracts the process of exporting Tukey
# results. This will be used in case I want to test for interactions/main
# effects not present in the first Tukey test function.


def exportTukeyResults(endFolder, parameter, indeVar, resultsTable):
    """
    Exports 2 files, a .txt file and an Excel file, that contains the raw Tukey
    results (containing difference in means, p-values, etc) from a Tukey test.
    These 2 files are copies of each other, with the results first exported
    into the .txt file and which will then be read and re-exported into an
    Excel file, which should be easier to handle and read

    Parameters
    ----------
    endFolder : str
        The name of the folder that these 2 files will end up in.
    parameter : str
        The parameter of interest. For enzymes, they will either be Vmax or Km.
        For litter chemistry, they will be the functional group of interest.
        This is the dependent variable that will be tested.
    indeVar : str
        The independent variable (or interaction) of interest. This will be
        used to delineate the different groups in this test.
    resultsTable : statsmodel table (format name unknown)
        The table containing Tukey results as outputted by statsmodel. Results
        will be read from this table and exported into the 2 export files.

    Returns
    -------
    None

    """
    folderPath = statsFolder/"Tukey posthoc"/endFolder
    txtName = "{0}, {1}.txt".format(parameter, indeVar)
    txtPath = folderPath/txtName
    if os.path.exists(txtPath) is False:
        # Creates and writes a text file of raw Tukey results if this file
        # doesn't exist yet
        txtFile = open(txtPath, "w")
        resultsRows = len(resultsTable)
        for i in range(resultsRows):
            row = resultsTable[i]
            for element in row:
                element = str(element)
                txtFile.write(element)
                txtFile.write(',')
            txtFile.write('\n')
        txtFile.close()
        # Creates an Excel file of raw Tukey results. Essentially
        # translating the .txt file over to the Excel file format
        resultsDF = pd.read_csv(txtPath, sep=",")
        excelName = "{0}, {1}.xlsx".format(parameter, indeVar)
        excelPath = folderPath/excelName
        if os.path.exists(excelPath) is False:
            resultsDF.to_excel(excelPath, index=False)
    return


# %%
# Conducting preliminary Tukey post-hoc on significant interactions/effects
# on litter chemistry. I may conduct more Tukey post-hoc on other effects later

# (1) Reading in the litter chemistry data and attaching it to the
# Michaelis-Menten parameters dataframe
litterChemPath = cwd/"Litter chemistry"/"Carbohydrates and Proteins FTIR.xlsx"
litterChem = pd.read_excel(litterChemPath)
litterChem.rename(columns={"functionalGroup": "Enzyme", "id": "ID"},
                  inplace=True)
litterChem["timePoint"] = litterChem["timePoint"].astype(str)
litterChem["Parameter"] = litterChem["Enzyme"]
parameters = pd.concat([parameters, litterChem])

# (2) Performing Tukey tests and exporting both raw Tukey test results and
# groups from each test result file

# Testing for carbohydrates
# Tukey("glycosidicBond", "glycosidicBond")
# groups("glycosidicBond")

# Tukey("C_O_stretching", "C_O_stretching")
# groups("C_O_stretching")

# Testing for total area of carbohydrate esters (carboEster1 + carboEster2.
# Pectins and hemicellulose contain esters and are subtypes of carbohydrates
# Tukey("carboEster", "carboEster")
# groups("carboEster")

# Testing for carboEster1 only
# Tukey("carboEster1", "carboEster1")

# Testing for carboEster2 only
# Tukey("carboEster2", "carboEster2")

# Testing for lipids
# Tukey("lipid", "lipid")
# groups("lipid")

# Testing for other alkanes
# Tukey("alkane", "alkane")
# groups("alkane")

# Testing for total amide area (amide1 + amide2). Amides are a functional
# group contained by proteins
# Tukey("amide", "amide")
# groups("amide")

# Testing for amide1 only
# Tukey("amide1", "amide1")

# Testing for amide2 only
# Tukey("amide2", "amide2")

# (3) Creating compact letter displays for the Tukey tests conducted so far on
# litter chemistry data
# for functionalGroup in functionalGroups:
#     litterChemTukeyFolder = statsFolder/"Tukey posthoc"/functionalGroup
#     files = os.listdir(litterChemTukeyFolder)

#     # Obtaining only the Excel raw Tukey results files, not the groups files
#     rawResultsFiles = []
#     for file in files:
#         if file.endswith(".xlsx") and "groups" not in file:
#             rawResultsFiles.append(file)

#     # Creating compact letter displays associated with each raw Tukey results
#     # file
#     for file in rawResultsFiles:
#         filePath = litterChemTukeyFolder/file
#         rawResults = pd.read_excel(filePath)
#         fileCLD = cld.main(rawResults)
#         # The cld.main() method returns a dataframe if the cld it creates,
#         # matches the raw Tukey results. If the cld doesn't match the raw
#         # Tukey results, then it doesn't return anything. This if statement
#         # is to check if the method returns a dataframe instead of nothing.
#         if type(fileCLD) is not None:
#             fileStripped = file.rstrip(".xlsx")
#             cldName = fileStripped + ", groups, annotated.xlsx"
#             cldPath = litterChemTukeyFolder/cldName
#             if os.path.exists(cldPath) is False:
#                 fileCLD.to_excel(cldPath, index=False)
#         elif type(fileCLD) is None:
#             print(file)
#             print("Failed to create CLD for these results")
#             print('\n')
"""Tukey's post-hoc rendered some main effects and interactions that were found
to be significant by ANOVAs to be insignificant. They are

Carbohydrate C-O stretching: time, precipitation
Carbohydrate total ester area: precipitation
carboEster1: time
carboEster2: time
Lipids: vegetation x precipitation
Alkanes: vegetation x precipitation
amide1: vegetation

As a result, I will conduct some additional Tukey's post-hoc on other main
effects/interactions that were still found to be significant by ANOVAs but
have not been tested yet by Tukey's post-hoc.
"""
# %%
# Purpose: conducting Tukey post-hoc on certain main effects/interactions on
# litter chemistry after preliminary Tukey's revealed that certain interactions
# are insignificant
# The main effects that will be tested for are
# Lipid: Vegetation

# (1) Conducting Tukey's post-hoc on Vegetation as a main effect on lipids
# and exporting them
# lipidDF = litterChem[litterChem.Enzyme == "lipid"]
# lipidResults = MC(lipidDF.value, lipidDF.Vegetation).tukeyhsd().summary()
# exportTukeyResults("lipid", "lipid", "Vegetation", lipidResults)

# (2) Creating the compact letter display for Vegetation as a main effect on
# lipids
# lipidTukeyFolder = statsFolder/"Tukey posthoc"/"lipid"
# lipidFiles = os.listdir(lipidTukeyFolder)
# for file in lipidFiles:
#     if file.endswith(".xlsx"):
#         # Obtaining only the Vegetation Excel Tukey results file
#         if "Vegetation" in file and "groups" not in file:
#             lipidTukeyPath = lipidTukeyFolder/file
#             lipidTukey = pd.read_excel(lipidTukeyPath)
#             lipidCLD = cld.main(lipidTukey)
#             if type(lipidCLD) is not None:
#                 fileStripped = file.rstrip(".xlsx")
#                 lipidCLDname = fileStripped + ", groups, annotated.xlsx"
#                 lipidCLDpath = lipidTukeyFolder/lipidCLDname
#                 lipidCLD.to_excel(lipidCLDpath, index=False)
#             else:
#                 print("CLD failed to be created for lipid, vegetation")
#                 break
"""And that's it. I can start making boxplots for litter chemistry data now"""
# %%
# Purpose: Performing Tukey's post-hoc on CAZyme domain relative abundance

# (1) Reading in the wrangled CAZyme data
CAZymePath = cwd/"CAZyme metagenomic data"/"Wrangled CAZyme gene domains.xlsx"
CAZ = pd.read_excel(CAZymePath)
CAZ.rename(columns={"Substrate": "Enzyme"}, inplace=True)
CAZ["timePoint"] = CAZ["timePoint"].astype(str)
parameters = pd.concat([parameters, CAZ])

# (2) Performing preliminary Tukey's post-hoc using the Tukey function I've
# written. Additional Tukey's post-hoc might be performed later. Also create
# compact letter displays for these Tukey results
# for substrate in substrates:
#     if substrate != "Total":
#         parameter = "Relative abundance"
#         folderName = substrate
#     elif substrate == "Total":
#         parameter = "Total CAZyme domains"
#         folderName = "Total CAZyme"
#     Tukey(substrate, parameter)  # Performing the Tukey test
#     # And now, creating the compact letter displays
#     folderPath = statsFolder/"Tukey posthoc"/folderName
#     folderContents = os.listdir(folderPath)
#     for content in folderContents:
#         # format of Tukey results in these folders are either a .txt file
#         # or an Excel .xlsx file
#         if content.endswith(".xlsx") and parameter in content:
#             resultsPath = folderPath/content
#             rawResults = pd.read_excel(resultsPath)
#             cldDF = cld.main(rawResults)
#             contentStrip = content.rstrip(".xlsx")
#             cldFileName = contentStrip + ", groups, annotated.xlsx"
#             cldPath = folderPath/cldFileName
#             if os.path.exists(cldPath) is False:
#                 cldDF.to_excel(cldPath, index=False)
"""I couldn't declare any interactions that were significant by ANOVA to be
insignificant by Tukey's. I'll make boxplots out of these results, and I'll
see if I still need to conduct more Tukey's"""
# %%
# Purpose: Perform additional Tukey's post-hoc (and make boxplots) on other
# main effects and interactions on CAZyme domain relative abundance.
"""The substrates that I will perform additional Tukey's post-hoc test on
includes

Cell wall: time, vegetation
Cellulose: time
Lignin: vegetation, precipitation
Polysaccharide: time, vegetation

I will also make boxplots using these additional Tukey results
"""

# (1) Performing additional Tukey's post-hoc
# Additional Tukey's on cell wall CAZyme domains
cellWallDF = parameters[parameters.Enzyme == "Cell_wall"]
cellWallTime = MC(cellWallDF.value, cellWallDF.timePoint).tukeyhsd().summary()
exportTukeyResults("Cell_wall", "Relative abundance", "timePoint",
                   cellWallTime)
cellWallVeg = MC(cellWallDF.value, cellWallDF.Vegetation).tukeyhsd().summary()
exportTukeyResults("Cell_wall", "Relative abundance", "Vegetation",
                   cellWallVeg)
# Additional Tukey's on cellulose CAZyme domains
celluloseDF = parameters[parameters.Enzyme == "Cellulose"]
celluloseTime = MC(celluloseDF.value,
                   celluloseDF.timePoint).tukeyhsd().summary()
exportTukeyResults("Cellulose", "Relative abundance", "timePoint",
                   celluloseTime)
# Additional Tukey's on lignin CAZyme domains
ligninDF = parameters[parameters.Enzyme == "Lignin"]
ligninVeg = MC(ligninDF.value, ligninDF.Vegetation).tukeyhsd().summary()
exportTukeyResults("Lignin", "Relative abundance", "Vegetation", ligninVeg)
ligninPpt = MC(ligninDF.value, ligninDF.Precip).tukeyhsd().summary()
exportTukeyResults("Lignin", "Relative abundance", "Precip", ligninPpt)
# Additional Tukey's on miscellaneous polysaccharide CAZyme domains
polysaccharideDF = parameters[parameters.Enzyme == "Polysaccharide"]
polysaccharideTime = MC(polysaccharideDF.value,
                        polysaccharideDF.timePoint).tukeyhsd().summary()
exportTukeyResults("Polysaccharide", "Relative abundance", "timePoint",
                   polysaccharideTime)
polysaccharideVeg = MC(polysaccharideDF.value,
                       polysaccharideDF.Vegetation).tukeyhsd().summary()
exportTukeyResults("Polysaccharide", "Relative abundance", "Vegetation",
                   polysaccharideVeg)

# (2) Creating compact letter displays for these additional Tukey's results
substratesRetest = ["Cell_wall", "Lignin", "Cellulose", "Polysaccharide"]
for substrate in substratesRetest:
    resultsFolder = statsFolder/"Tukey posthoc"/substrate
    allFiles = os.listdir(resultsFolder)
    rawResults = []
    cldFiles = []

    # Obtaining 2 lists, 1 of only raw results, 1 of only the compact letter
    # displays created from those results
    for file in allFiles:
        strippedFile = file.rstrip(".xlsx")
        if file.endswith(".xlsx") and "annotated" not in file:
            rawResults.append(strippedFile)
        elif "annotated" in file:
            cldFiles.append(strippedFile)
    """Now, starting a process in which I obtain only the raw results from
    these additional Tukey post-hoc tests"""

    # Obtaining a list of raw results from prior Tukey's tests, not these new
    # Tukey's tests. I will then find the difference between this list
    # and the list of all raw results
    priorResults = []
    for cldFile in cldFiles:
        for rawFile in rawResults:
            if rawFile in cldFile:
                priorResults.append(rawFile)

    # Obtaining the set of the results from the additional Tukey's post-hoc
    # using a set operation
    rawResults = set(rawResults)
    priorResults = set(priorResults)
    additionalResults = rawResults - priorResults  # Set operation in which
    # you find the items that are only unique to the first set, rawResults,
    # so files that are in rawResults but not in priorResults, which would
    # be the newly created Tukey results files

    # Creating and exporting a compact letter display of each additional
    # results file
    for file in additionalResults:
        cldName = file + ", groups, annotated.xlsx"
        resultsExcel = file + ".xlsx"
        resultsPath = resultsFolder/resultsExcel
        results = pd.read_excel(resultsPath, na_filter=False)
        cldResults = cld.main(results)
        cldPath = resultsFolder/cldName
        if os.path.exists(cldPath) is False:
            cldResults.to_excel(cldPath, index=False)
# %%
print(dt.now() - start)
