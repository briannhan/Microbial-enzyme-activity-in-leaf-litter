# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 23:24:06 2020

@author: Brian Chung
The purpose of this script is to calculate the enzyme activity of the litter
from Loma Ridge. This script processes enzyme activity & litter dry weight
data.

I delineate this script into sections using the following characters: # %%
These sections are visualized in the Spyder IDE of the Anaconda distribution
of Python. Each section has a particular purpose, and several tasks are
numbered to accomplish that purpose. The number in the tasks doesn't reflect
the order at which they are done. For examples, if I list tasks (1), (2), and
(3), I may complete task (3) before task (2). However, all the tasks in each
section will be completed to carry out the purpose of that section.
"""
# %%
# Importing necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as py
import os
from pathlib import Path

# %%
# Making file & folder paths to read in data
workingDirectory = Path(os.getcwd())
enzymeFolder = "Enzyme activity data"
dryMassFolder = "Litter dry weights"
enzymeFolderPath = workingDirectory/enzymeFolder
enzymeFiles = os.listdir(enzymeFolderPath)
dryMassFolderPath = workingDirectory/dryMassFolder
dryMassFiles = os.listdir(dryMassFolderPath)

# %%
# Processing dry weight data before using it for enzyme activity calculations
# Processing includes:
# (1) renaming columns
# (2) Filling in time points;
# (3) Adding columns for wet and dry assay mass and proportion of leaf litter
# as dry mass and calculating proportion dry assay mass; and
# (4) isolating timepoints T0, T3, and T5 (T6 will be isolated from a different
# spreadsheet)

# The following line makes pandas print all columns
pd.set_option("display.max_columns", None)

# Reading in dry mass file
dryMassFilePath = dryMassFolderPath/dryMassFiles[0]
dryDFfull = pd.read_excel(dryMassFilePath)

# (1) Renaming columns
oldCols = dryDFfull.columns.tolist()
newCols = ["Envelope mass (g)", "Envelope + wet (g)", "Envelope + dry (g)",
           "Wet litter mass (g)", "Dry litter mass (g)"]

newColumnsDict = {oldCols[n + 2]: newCols[n] for n in range(len(newCols))}
dryDFfull = dryDFfull.rename(columns=newColumnsDict)


# (2) Filling the time column with time points
for index, row in dryDFfull.iterrows():
    currentTimepointType = type(row["Time"])
    if currentTimepointType != str:
        previousRow = dryDFfull.iloc[index - 1]
        dryDFfull.loc[index, "Time"] = previousRow["Time"]

# (3) Adding columns & calculating dry assay mass
dryDFfull["Wet assay (g)"] = 0.4
dryDFfull["Dry proportion (%)"] = 100*(dryDFfull["Dry litter mass (g)"]
                                       / dryDFfull["Wet litter mass (g)"])
dryDFfull["Dry assay (g)"] = (dryDFfull["Wet assay (g)"]
                              * dryDFfull["Dry proportion (%)"]/100)

# (4) Isolating timepoints T0, T3, T5
timepoints = dryDFfull.groupby("Time")["Time"].count().index.tolist()
timepoints = [timepoints[0], timepoints[3], timepoints[5]]
T0_3_5Bool = dryDFfull["Time"].isin(timepoints)
dryDFfull = dryDFfull[T0_3_5Bool]

oldCols = dryDFfull.columns.tolist()
colsToDrop = oldCols[2:-1]
dryDFprocessed = dryDFfull.drop(labels=colsToDrop, axis=1)
# %%
# I'm going to test with a single enzyme file first as a proof of concept. In
# other words, I'm going to do some pre-processing for a single enzyme file
# in this section before going on to calculate enzyme activities in the
# following sections.
# The tasks for this section are:
# (1) Process long sample names in the enzyme file;
# (2) Check to see if the names in the enzyme files match the names in the dry
# weights files; if there isn't a match, then some samples are missing or some
# samples are redone
# (3) rearrange control readings (substrate, homogenate, and quench controls)
# to ease calculations
# (4) add additional data to the enzyme data frame, which consists of dry assay
# mass and vegetation & precipitation treatments

# Obtaining a list of samples to check if plate data contain samples
samples = dryDFfull.groupby("ID")["ID"].count().index.tolist()

# Reading in T0 black plates
ogEnzymeCols = ["Well", "Long sample name", "Fluorescence"]
T0BlackPath = enzymeFolderPath/enzymeFiles[0]
T0Black = pd.read_csv(T0BlackPath, '\t', header=None, names=ogEnzymeCols)

# (1) Processing long sample names. The long sample names have a lot of junk in
# them so I'm processing them to keep (1) the assay date and (2) information
# for whether the plate is a buffer of a sample plate. If the plate is a sample
# plate, then it would contain the sample ID
T0Black["Long sample name"] = T0Black["Long sample name"].str.split("_")
for index, row in T0Black.iterrows():
    longSampleName = row["Long sample name"]
    T0Black.loc[index, "Assay date"] = longSampleName[0]
    if "B" in longSampleName:
        T0Black.loc[index, "ID"] = "B"
    elif "X" in longSampleName[2]:
        T0Black.loc[index, "ID"] = longSampleName[2]
T0Black = T0Black.drop(labels="Long sample name", axis=1)

# (2) Checking to see if enzyme activity data contains all samples
T0BlackCounts = T0Black.groupby("ID")["ID"].count()/96
# If any of the counts of the samples return 2 or greater, then that means
# that that sample was assayed at least twice. But what if one sample has
# a different name than what the name should be?
T0BlackSamples = T0BlackCounts.index.tolist()
samplesNotInT0Black = []
for sample in samples:
    if sample not in T0BlackSamples:
        samplesNotInT0Black.append(sample)
# if samplesNotInT0Black is empty, that means that all the sample names in
# the dry weights data frame match the sample names in the enzyme data, and so
# no samples are mis-named or missing


# (4) Adding additional data to the enzyme data frame
# Adding dry assay mass
T0Bool = dryDFprocessed["Time"].isin(["T0 (November 30, 2017)"])
dryDF_T0 = dryDFprocessed[T0Bool]
T0Black = pd.merge(T0Black, dryDF_T0, how="left", on="ID")
T0Black = T0Black.drop(labels="Time", axis=1)
# be sure to read up on pandas documentation to figure out how the merge
# method works. Go to section 2.7.2 of the guide on pandas version 1.1.4, page
# 448 (462 on the online pdf reader)
# Adding treatments
for index, row in T0Black.iterrows():
    rowID = row["ID"]
    if rowID != "B":
        precip = rowID[-2]
        plot = int(rowID[:-3])
        if plot <= 24:
            T0Black.loc[index, "Vegetation"] = "Grassland"
        elif plot > 24:
            T0Black.loc[index, "Vegetation"] = "CSS"
        if precip == "X":
            T0Black.loc[index, "Precip"] = "Ambient"
        elif precip == "R":
            T0Black.loc[index, "Precip"] = "Reduced"

# (3) Rearranging control readings. I intend that, for each sample, the
# absorbance in a particular well will also correspond with the absorbance
# reading in the buffer plate at that particular well. There will be 2 columns
# one column to hold absorbance for a particular well of the sample while the
# second column holds the absorbance for the same well for the buffer plate.
# I will also manipulate the quench & homogenate control readings to make sure
# that they are side by side with the sample readings

bufferDF = T0Black[T0Black["ID"] == "B"]
T0Black = T0Black[T0Black["ID"] != "B"]
oldCols = T0Black.columns.tolist()
colsToDrop = oldCols[-4:]
bufferDF = bufferDF.drop(labels=colsToDrop, axis=1)
bufferDF = bufferDF.rename(mapper={"Fluorescence": "BufferReading"}, axis=1)

T0Black = pd.merge(T0Black, bufferDF, how="inner", on=["Well", "Assay date"])
# sortCols = ["ID"]
# T0Black = T0Black.sort_values(by=sortCols, axis=1)
