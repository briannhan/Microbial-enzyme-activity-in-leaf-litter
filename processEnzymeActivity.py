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
oxiErrorLabels = ["PPO rep 1", "PPO rep 2", "PER rep 1", "PER rep 2"]
T0hydroProcessing = T0processing.drop(labels=oxiErrorLabels, axis=1)
hydroEnzymes = T0hydroProcessing.columns.tolist()[1:]
initialShape = T0hydroData.shape
print("Initial shape is", initialShape)
for index, row in T0hydroProcessing.iterrows():
    for enzyme in hydroEnzymes:
        if row[enzyme] == "o":  # substrate inhibition
            sampleDF = T0hydroData[T0hydroData["ID"] == row["ID"]]
            dfToProcess = sampleDF[sampleDF["Enzyme"] == enzyme]
            indexMax = dfToProcess["Activity"].idxmax()
            dfIndices = dfToProcess.index.tolist()
            indicesToDrop = np.arange(dfIndices[0], indexMax + 1)
            # print(row["ID"], enzyme, "all indices:", dfIndices, "max index:",
            #       indexMax, "indices to drop:", indicesToDrop)
            T0hydroData = T0hydroData.drop(index=indicesToDrop)
finalShape = T0hydroData.shape
print("Final shape is", finalShape)
# %%
# Purpose: Fit the single-substrate, single-enzyme approximation of equilibrium
# chemistry (ECA) by Tang and Riley 2013 to the processed T0 hydrolase activity
# data and store the parameters in a new dataframe


def ECA(k2, E, Km):
    """
    This is the single-substrate, single-enzyme equilibrium chemistry
    approximation (ECA) developed by Tang and Riley 2013. The units of the
    input data have been normalized accordingly to ensure that the units of
    the parameters in this function are also normalized properly when this
    function is fitted to the data.

    Parameters
    ----------
    k2 : TYPE
        The rate constant of the 2nd step in equilibrium chemistry. Units are
        micromolar g^-1 s^-1.
    E : TYPE
        The concentration of enzymes in a well. Units are micromolar g^-1.
    Km : TYPE
        The concentration of substrates in a well. Units are micromolar g^-1.

    Returns
    -------
    Fitted enzyme activity
    """
    V = 
# %%
print(datetime.now() - startTime)
