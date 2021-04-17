# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:18:28 2021

@author: Brian Chung
The purpose of this script is to process enzyme activity that had been
calculated using the calculateEnzymeActivity.py script. This script uses an
Excel spreadsheet file I've made containing my labels for how I processed the
errors I saw when looking at the unprocessed enzyme activity graphs.

This script processes enzyme activity by time point.
After each processing step, the data will be fitted to the single-substrate,
single-enzyme formulation of an approximation of equilibrium chemistry (ECA)
according to Tang and Riley 2013. The fitting of the data will be constrained
so that values of k2 (the maximum product genesis rate), Km (Michaelis-Menten
constant), and enzyme concentration are all positive.
"""
import os
import pandas as pd
from pathlib import Path
from datetime import datetime
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as py
import enzymeWrangling as ew
startTime = datetime.now()

# Accessing the necessary files, which are (1) Excel file of unprocessed,
# calculated enzyme activity and (2) Excel file that details the errors &
# modifications I should make of the unprocessed enzyme activity data.
cwd = Path(os.getcwd())
cwdDirs = os.listdir(cwd)
unprocessedGraphsPath = cwd/"Unprocessed activity graphs"
unprocessedGraphsContents = os.listdir(unprocessedGraphsPath)
processingTypesPath = unprocessedGraphsPath/"Samples with errors.xlsx"
''' "types" here is short for "types of processing" that I'll do, which I've
# annotated in the "Samples with errors.xlsx" Excel file'''
processingTypes = pd.ExcelFile(processingTypesPath)
processingSheets = processingTypes.sheet_names

enzymeActivityFolder = cwd/"Enzyme activity data"
enzymeActivityFiles = os.listdir(enzymeActivityFolder)
activityDataPath = enzymeActivityFolder/"Unprocessed Enzyme Activity.xlsx"
activityData = pd.ExcelFile(activityDataPath)
activitySheets = activityData.sheet_names
# I'm gonna process T0 first as a test, and then write methods to generalize
# and abstract the processing steps in enzymeWrangling.py
# %%
# Purpose: Processing T0 Hydrolase activity data
# (1) Reading in files concerning T0 Hydrolases from the Excel file where I
# annotated the processing necessary and from the Excel file containing the
# calculated enzyme activit data
# (2) Set all negative values to equal 0
# (3) Process substrate inhibition by removing data points with lower activity
# but higher substrate concentration than the data point with the highest
# enzyme activity

# (1) Reading in files concerning T0 Hydrolases from the Excel file where I
# annotated the processing necessary and from the Excel file containing the
# calculated enzyme activit data
T0hydroData = pd.read_excel(io=activityData, sheet_name=activitySheets[0])
T0processing = pd.read_excel(processingTypes, sheet_name=processingSheets[0])
T0processing = T0processing.drop(labels="Unnamed: 12", axis=1)

# (2) Set all negative values to equal 0
T0hydroCols = T0hydroData.columns.tolist()
T0hydroProcessInd = T0hydroData.index[T0hydroData["Activity"] < 0].tolist()
T0hydroData.loc[T0hydroProcessInd, "Activity"] = 0

# (3) Process substrate inhibition by removing data points with lower activity
# but higher substrate concentration than the data point with the highest
# enzyme activity
samples = T0processing["ID"].tolist()
oxiErrorLabels = ["PPO rep 1", "PPO rep 2",
                  "PER rep 1", "PER rep 2",]
T0hydroProcessing = T0processing.drop(labels=oxiErrorLabels, axis=1)
hydroEnzymes = T0hydroProcessing.columns.tolist()[1:]
initialShape = T0hydroData.shape
print(T0hydroData)
print("Initial shape is", initialShape)
for index, row in T0hydroProcessing.iterrows():
    for enzyme in hydroEnzymes:
        if row[enzyme] == "o":  # substrate inhibition
            sampleDF = T0hydroData[T0hydroData["ID"] == row["ID"]]
            dfToProcess = sampleDF[sampleDF["Enzyme"] == enzyme]
            indexMax = dfToProcess["Activity"].idxmax()
            dfIndices = dfToProcess.index.tolist()
            indicesToDrop = np.arange(dfIndices[0], indexMax)
            # print(row["ID"], enzyme, "all indices:", dfIndices, "max index:",
            #       indexMax, "indices to drop:", indicesToDrop)
            T0hydroData = T0hydroData.drop(index=indicesToDrop)
finalShape = T0hydroData.shape
print("Final shape is", finalShape)
# %%
# Purpose: Fit the single-substrate, single-enzyme approximation of equilibrium
# chemistry (ECA) by Tang and Riley 2013 to the processed T0 hydrolase activity
# data and store the parameters in a new dataframe
# Tasks:
# (1) Create ECA function
# (2) Fit function to T0 hydrolase dataframe

# (1) Create ECA function


def ECA(S, k2, E, Km):
    """
    This is the single-substrate, single-enzyme equilibrium chemistry
    approximation (ECA) developed by Tang and Riley 2013. The units of the
    input data have been normalized accordingly to ensure that the units of
    the parameters in this function are also normalized properly when this
    function is fitted to the data.

    Parameters
    ----------
    S : TYPE
        The concentration of substrates in a well. Units are micromolar g^-1.
    k2 : TYPE
        The rate constant of the 2nd step in equilibrium chemistry. Units are
        L hr^-1.
    E : TYPE
        The concentration of enzymes in a well. Units are micromolar g^-1.
    Km : TYPE
        Michaelis-Menten constant. The equilibrium constant that is the ratio
        between the rate of breakdown of the enzyme-substrate complex and the
        formation of the enzyme-substrate complex. Units are micromolar g^-1.

    Returns
    -------
    Fitted enzyme activity. Units are micromole hr^-1
    """
    return (k2*E*S)/(Km + E + S)


# (2) Fit function to T0 hydrolase dataframe
samplesDict = {"ID": samples}
T0k2 = pd.DataFrame(data=samplesDict)
T0enzymeConcen = pd.DataFrame(data=samplesDict)
T0Km = pd.DataFrame(data=samplesDict)
test4LXX = T0hydroData[T0hydroData["ID"] == "4LXX"]
test4LXX_AG = test4LXX[test4LXX["Enzyme"] == "AG"]
params1, paramCov1 = curve_fit(ECA, test4LXX_AG["NormSubConcen"],
                               test4LXX_AG["Activity"],
                               bounds=(0, np.inf))
params2, paramCov2 = curve_fit(ECA, test4LXX_AG["NormSubConcen"],
                               test4LXX_AG["Activity"])
# testFig = py.figure(figsize=(8, 6))
# subplot1 = testFig.add_subplot(1, 1, 1)
# py.title("4LXX AG test")
# py.xlabel("Normalized substrate concentration (micromolar/g)")
# py.ylabel("Normalized enzyme activity (micromole/hr)")
# py.plot("NormSubConcen", "Activity", data=test4LXX_AG, marker="o",
#         linestyle="-", color="b",)
# testXvals = np.linspace(0, test4LXX_AG["NormSubConcen"].max())
# testY1 = ECA(testXvals, params1[0], params1[1], params1[2])
# testY2 = ECA(testXvals, params2[0], params2[1], params2[2])
# py.plot(testXvals, testY1, "-g")
# py.plot(testXvals, testY2, "-r")
# T0k2.loc["4LXX", "AG"] = params[0]
# T0enzymeConcen.loc["4LXX", "AG"] = params[1]
# T0Km.loc["4LXX", "AG"] = params[2]
'''So, I noticed something weird. When I didn't set any boundaries, I get a
specific set of parameter values. When I set the boundaries of the parameter
values so that they are between 0 and positive infinity, I get a different set
of parameter values even though the first set of values are all positive. Why?

And also, when I plotted both sets of parameters, both of them fit the data
equally well to each other, so much so that the values they predict are
indistinguishable from each other unless you zoom in very, very closely.
So, shit. Wonder what I should do now.

Still, setting the bounds of parameter values to be between 0 & positive
infinity is, I think, very important, so I'll keep doing that.
'''

# for sample in samples:
#     sampleDF = T0hydroData[T0hydroData["ID"] == sample]
#     sampleIndex = T0k2[T0k2["ID"] == sample].index
#     for enzyme in hydroEnzymes:
#         enzymeDF = sampleDF[sampleDF["Enzyme"] == enzyme]
#         try:
#             params, paramCov = curve_fit(ECA, enzymeDF["NormSubConcen"],
#                                          enzymeDF["Activity"])
#             T0k2.loc[sampleIndex, enzyme] = params[0]
#             T0enzymeConcen.loc[sampleIndex, enzyme] = params[1]
#             T0Km.loc[sampleIndex, enzyme] = params[2]
#         except RuntimeError:
#             print("Can't fit sample {0:}, enzyme {1:}".format(sample,
#                                                               enzyme))
#             T0k2.loc[sampleIndex, enzyme] = "can't fit"
#             T0enzymeConcen.loc[sampleIndex, enzyme] = "can't fit"
#             T0Km.loc[sampleIndex, enzyme] = "can't fit"
# T0k2Melt = pd.melt(T0k2, "ID", var_name="Enzyme", value_name="k2")

"""This nested for loop failed at the AG enzyme of sample 25LRX before I added
in the try/accept statement. Might need to throw that sample out. Also, Steve
suggested that I use Michaelis-Menten instead so that I'd be less likely to run
into this identifiability problem and also to make it easier to compare to
prior work.

Using Michaelis-Menten means that I'm gonna have to rework my hypotheses,
though. Boo.

I've went ahead and copied and pasted the code above into enzymeWrangling.py
and reworking this code in order to apply it to both hydrolytic and oxidative
enzymes and to oxidative enzymes in T5, in which sample 47 was assayed twice.
"""
# %%
# Purpose: Obtain Michaelis-Menten parameters for T0 hydrolase data
T0params = ew.nonlinRegress(T0hydroData, "H")
# %%
print(datetime.now() - startTime)
