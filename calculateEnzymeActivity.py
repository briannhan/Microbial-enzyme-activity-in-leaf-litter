# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 23:24:06 2020

@author: Brian Chung
The purpose of this script is to calculate the enzyme activity of the litter
from Loma Ridge. This script processes enzyme activity & litter dry weight
data. Refer to the 2020-11-10-enzyme-descriptive-metadata.docx Word Document
for a description of how assays are carried out.

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
import enzymeWrangling as ew
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

# Purpose: begin pre-processing enzyme data from T0 black
# The tasks for this section are:
# (1) Process long sample names in the enzyme file;
# (2) Check to see if the names in the enzyme files match the names in the dry
# weights files; if there isn't a match, then some samples are missing or some
# samples are redone
# (3) add additional data to the enzyme data frame, which consists of dry assay
# mass and vegetation & precipitation treatments

# Obtaining a list of samples to check if plate data contain samples
samples = dryDFfull.groupby("ID")["ID"].count().index.tolist()

# Reading in T0 black plates
ogEnzymeCols = ["Well", "Long sample name", "Assay"]
T0BlackPath = enzymeFolderPath/enzymeFiles[0]
T0Black = pd.read_csv(T0BlackPath, '\t', header=None, names=ogEnzymeCols)

# (1) Processing long sample names. The long sample names have a lot of junk in
# them so I'm processing them to keep (1) the assay date and (2) information
# for whether the plate is a buffer of a sample plate. If the plate is a sample
# plate, then it would contain the sample ID
# I'm also splitting wells into letters that represent rows and numbers that
# represent columns on a black plate
T0Black["Long sample name"] = T0Black["Long sample name"].str.split("_")
for index, row in T0Black.iterrows():
    longSampleName = row["Long sample name"]
    T0Black.loc[index, "Assay date"] = longSampleName[0]
    if "B" in longSampleName:
        T0Black.loc[index, "ID"] = "B"
    elif "X" in longSampleName[2]:
        T0Black.loc[index, "ID"] = longSampleName[2]
    well = row["Well"]
    wellRow = well[0]
    wellColumn = int(well[1:])
    T0Black.loc[index, "PlateRow"] = wellRow
    T0Black.loc[index, "PlateCol"] = wellColumn
T0Black = T0Black.drop(labels="Long sample name", axis=1)
T0Black = T0Black[T0Black["PlateCol"] <= 10]

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


# (3) Adding additional data to the enzyme data frame
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
        T0Black.loc[index, "Plot"] = plot
# %%
# Purpose: finish pre-processing of T0 Black by rearranging control readings
# I intend that, for each sample, the
# fluorescence in a particular well will also correspond with the fluorescence
# reading in the buffer plate at that particular well. There will be 2 columns
# one column to hold fluorescence for a particular well of the sample while the
# second column holds the fluorescence for the same well for the buffer plate.
# I will also manipulate the quench & homogenate control readings to make sure
# that they are side by side with the sample readings
# Tasks of this section:
# (1) Manipulating substrate control
# (2) Manipulating quench control & standard fluorescence
# (3) Manipulating homogenate control
# Note on running this cell: You cannot run this cell all by itself. You must
# run the prior cell (which reads in T0Black) before you can run this cell.

# (1) Manipulating substrate control data by placing the buffer & sample
# readings side by side. The buffer readings for columns 1-7 on the black plate
# will serve as the substrate control.
bufferDF = T0Black[T0Black["ID"] == "B"]
T0Black = T0Black[T0Black["ID"] != "B"]
oldCols = T0Black.columns.tolist()
colsToDrop = ["ID", "Dry assay (g)", "Vegetation", "Precip", "Plot"]
bufferDF = bufferDF.drop(labels=colsToDrop, axis=1)
bufferDF = bufferDF.rename(mapper={"Assay": "BufferReading"}, axis=1)
subCtrlMergeLabels = ["Well", "Assay date", "PlateRow", "PlateCol"]
T0Black = pd.merge(T0Black, bufferDF, how="inner", on=subCtrlMergeLabels)
T0Black = T0Black.sort_values(by=["PlateCol", "Plot"])

# Calculations of enzyme activity follow from the German et al 2011 paper:
# Optimization of hydrolytic and oxidative enzyme methods for ecosystem studies
# Refer to this paper for reference if necessary

# (2) Manipulating standard fluorescence and quench control readings.
# columns 8 & 9 of the BUFFER plate represent standard fluorescence, while the
# same columns in a SAMPLE plate represent quench fluorescence
AMC_DF = T0Black[T0Black["PlateCol"] == 8]
AMC_DF = AMC_DF.drop(labels="Well", axis=1)
MUB_DF = T0Black[T0Black["PlateCol"] == 9]
MUB_DF = MUB_DF.drop(labels="Well", axis=1)
homCtrlDF = T0Black[T0Black["PlateCol"] == 10]
T0Black = T0Black[T0Black["PlateCol"] <= 7]

# Setting black plate columns to their respective standards.
MUB_DF1 = MUB_DF.copy()
MUB_DF2 = MUB_DF.copy()
MUB_DF3 = MUB_DF.copy()
MUB_DF4 = MUB_DF.copy()
MUB_DF5 = MUB_DF.copy()
MUB_DF7 = MUB_DF.copy()

standardFrames = [MUB_DF1, MUB_DF2, MUB_DF3, MUB_DF4, MUB_DF5, AMC_DF, MUB_DF7]
for n in range(len(standardFrames)):
    currentDF = standardFrames[n]
    currentDF["PlateCol"] = n + 1

# Merging AMC_DF & MUB_DF into a single dataframe of standard fluorescence
# & quench control. This subsequent dataframe will be merged back into
# T0Black
standardDF = pd.concat(objs=standardFrames, axis=0)
standardDF = standardDF.sort_values(by=["PlateCol", "Plot"])
# This is the right way to sort any dataframe that's derived from the enzyme
# data to ensure that the sorted data frame looks like the original file

# Renaming columns in the standard data frame and dropping redundant columns.
# I deem certain columns as redundant because this dataframe will be merged
# back into T0 Black, which already contains the information in the redundant
# columns
newColumnsDict = {"Assay": "QuenchCtrl",
                  "BufferReading": "StanFluo"}
stanColsToDrop = ["Dry assay (g)", "Vegetation", "Precip"]
standardDF = standardDF.rename(columns=newColumnsDict)
standardDF = standardDF.drop(labels=stanColsToDrop, axis=1)

# Merging standard dataframe back into T0Black
stanMergeLabels = ["Assay date", "ID", "PlateRow", "PlateCol", "Plot"]
T0Black = pd.merge(T0Black, standardDF, how="inner", on=stanMergeLabels)
T0Black = T0Black.sort_values(by=["Plot", "PlateCol"])  # sorting T0Black
T0Black = T0Black.rename(columns={"BufferReading": "SubCtrl"})


""" The following memo is for me to follow, so it might not make sense to
random online viewers of this project. This memo mentions a script by my PI
that I will hereafter refer to as "Steve's script". Random online viewers of
this project will, of course, not know what this script is. To summarize, the
script is an example I follow to develop this calculateEnzymeActivity.py
script. Thus far, I've been wondering whether I should use column 10 of both
the sample plates and buffer plates, but it looks like my PI only uses column
10 of the sample plate, so I will also use only column 10 of the sample plate.
Column 10 of the buffer plates seem useless, as their fluorescence tend to be
0 are close to 0 and so does not change enzyme activities much when added or
subtracted.

So I read Steve's script. Line 55 is where he calculates the enzyme activity
for black plates. He only uses column 10 of the sample plate, not of the
buffer plate, if I'm interpreting his script correctly. He classifies
black plate columns using the EnzymePlateInfo.xlsx spreadsheet. Column 10 of
the buffer plate is classified as "Buffer" while for a sample plate it's
classified as "HomBlank". While "Buffer" was used in line 32 of the script,
it wasn't actually used in line 55. On the other hand, HomBlank is used in
line 55, so this makes me think that Steve only uses column 10 of the sample
plate

So column 10 of only the SAMPLE plates will serve as the homogenate control,
controlling for any fluorescence that the homogenate itself might have"""
# (3) Manipulating homogenate control data. I'll also be dropping redundant
# columns
homCtrlDF = homCtrlDF.rename(mapper={"Assay": "HomCtrl"}, axis=1)
homCtrlColsToDrop = ["Well", "Assay date", "PlateCol", "Dry assay (g)",
                     "Vegetation", "Precip", "BufferReading"]
homCtrlDF = homCtrlDF.drop(labels=homCtrlColsToDrop, axis=1)
homCtrlMergeLabels = ["ID", "Plot", "PlateRow"]
T0Black = pd.merge(T0Black, homCtrlDF, how="inner", on=homCtrlMergeLabels)
T0Black = T0Black.sort_values(by=["Plot", "PlateCol"])

"""Let's talk about how to plot the data. For each sample, let's make a figure
with 7 subplots, one for each enzyme. Let's create a function where each call
plots a single sample. To plot all samples, I'll use a for loop to loop through
all the sample names and for each sample name, I'll plot that sample.

And now, the body of the function. The function will filter the T0Black data
frame down to a single sample. The function will also create 7 data frames,
1 for each enzyme. It'll then call the plot() and scatter() methods of
matplotlib.pyplot 7 times, 1 for each enzyme, where the x-values would be the
substrate concentrations and the y-values would be the calculted enzyme
activities"""
# %%
# Purpose: Calculate hydrolytic enzyme activity of T0Black using formulas that
# convert fluorescence to enzyme activity in German et al 2011
# Tasks:
# (1) Making a dataframe of enzyme plate information, containing enzyme name,
# standard for the enzyme, substrate concentrations, and standard amounts
# (2) Calculating quench coefficients
# (3) Calculating emission coefficients
# (4) Calculating net fluorescence
# (5) Calculating enzyme activity

# (1) Making dataframe to hold enzyme name, standard for the enzyme, substrate
# concentrations, and standard amounts. This dataframe essentially holds
# certain plate information and is called plateInfo.
# plateInfo will, at first, hold the enzyme name and amounts of standards in
# columns 8 & 9. I'm creating another dataframe of enzyme
# concentrations to merge back into plateInfo.
AMCamount = 62.5*125/1000
"""Amount of MUB standard in standard and quench control wells.
62.5 micromolar x 125 microliter / 1000 (conversion factor)
units: nanomoles"""
MUBamount = 25*125/1000
"""Amount of AMC standard in standard and quench control wells.
25 micromolar x 125 microliter / 1000 (conversion factor)
units: nanomoles"""
plateCols = np.linspace(start=1, stop=7, num=7).tolist()
enzymeName = ["AG", "AP", "BG", "BX", "CBH", "LAP", "NAG"]
standardAmount = [MUBamount, MUBamount, MUBamount, MUBamount, MUBamount,
                  AMCamount, MUBamount]
plateInfoDic = {"PlateCol": plateCols, "Enzyme": enzymeName,
                "StanAmt": standardAmount}
plateInfo = pd.DataFrame(plateInfoDic)
# Making dataframe of hydrolytic enzyme concentrations to merge back into
# plateInfo.
AGnames = 8*["AG"]
APnames = 8*["AP"]
BGnames = 8*["BG"]
BXnames = 8*["BX"]
CBHnames = 8*["CBH"]
LAPnames = 8*["LAP"]
NAGnames = 8*["NAG"]
enzymesLists = [AGnames, APnames, BGnames, BXnames,
                CBHnames, LAPnames, NAGnames]
longEnzymesNames = [name for outerList in enzymesLists for name in outerList]
subProps = np.geomspace(1, 1/128, num=8)
# Units of substrate concentrations are in micromolar
AGconcen = (500*subProps).tolist()
APconcen = (2000*subProps).tolist()
BGconcen = (500*subProps).tolist()
BXconcen = (500*subProps).tolist()
CBHconcen = (250*subProps).tolist()
LAPconcen = (500*subProps).tolist()
NAGconcen = (1000*subProps).tolist()
subConcenLists = [AGconcen, APconcen, BGconcen, BXconcen, CBHconcen,
                  LAPconcen, NAGconcen]
subConcen = [concen for outer in subConcenLists for concen in outer]
plateRow = 7*list("ABCDEFGH")
blackSubConcenDict = {"PlateRow": plateRow, "Enzyme": longEnzymesNames,
                      "SubConcen": subConcen}
blackSubConcenDF = pd.DataFrame(blackSubConcenDict)
plateInfo = pd.merge(plateInfo, blackSubConcenDF, how="inner", on="Enzyme")
T0Black = pd.merge(T0Black, plateInfo, on=["PlateCol", "PlateRow"])


# (2) Calculating Quench Coefficients. Each quench coefficient is specific
# to a particular row of a sample plate. Quench Coefficients are unitless
T0Black["QuenchCoef"] = ((T0Black["QuenchCtrl"] - T0Black["HomCtrl"])
                         / T0Black["StanFluo"])

# (3) Calculating Emission Coefficients. Each coefficient is specific to a row
# of a plate. Units: nmol^-1
assayVol = 0.250  # mL. Volume of assay well, which consists of
# substrate (0.125 mL) + homogenate (0.125 mL)
stanVol = 0.250  # mL. Volume of quench and standard wells, which
# consists of standard (0.125 mL) +
# either homogenate(sample plate, 0.125 mL) or buffer(buffer plate, 0.125 mL)
T0Black["EmisCoef"] = ((T0Black["StanFluo"]*stanVol)
                       / (T0Black["StanAmt"]*assayVol))

# (4) Calculating net fluorescence. Units are fluorescence units, which is
# essentially meaningless and so can be considered as unitless
T0Black["NetFluo"] = (((T0Black["Assay"] - T0Black["HomCtrl"])
                       / T0Black["QuenchCoef"]) - T0Black["SubCtrl"])

# (5) Calculating enzyme activity. Initial units are nmol L^-1 g^-1 h^-1. Final
# units will be micromole L^-1 g^-1 h^-1. This enzyme activity had been
# normalized for the mass of litter used in the assay
incubaTime = 4  # Hours. Each plate was left to sit for 4 hours
bufferVol = 150  # mL. Volume of buffer used in preparing a single homogenate
homVol = 0.125  # mL. Volume of homogenate that is pipetted into each well
# in columns 1-7 in each sample plate.
T0Black["Activity"] = ((T0Black["NetFluo"]*bufferVol)
                       / (T0Black["EmisCoef"]*homVol*incubaTime
                          * T0Black["Dry assay (g)"]))
T0Black["Activity"] = T0Black["Activity"]/1000
T0Black = T0Black.sort_values(by=["Plot", "PlateCol"])
# %%
# Purpose: plotting hydrolytic enzyme activity of T0. I won't be doing any
# nonlinear regression to fit the ECA model to the enzyme activity just yet
py.style.use("dark_background")
T0graphsFolder = (workingDirectory/"Unprocessed activity graphs"
                  / "T0"/"Hydrolytic enzymes")
for sample in samples:
    sampleDF = T0Black[T0Black["ID"] == sample]
    vegetation = sampleDF.groupby("Vegetation")["Vegetation"].count().index
    vegetation = vegetation.tolist()
    vegetation = vegetation[0]
    precip = sampleDF.groupby("Precip")["Precip"].count().index.tolist()
    precip = precip[0]
    figTitle = "{0:}, {1:}, {2:}, T0 Black".format(sample, vegetation, precip)
    py.figure(num=figTitle, figsize=(25, 10))
    for enzyme in enzymeName:
        if enzyme == "NAG":
            plotIndex = 8
        else:
            plotIndex = enzymeName.index(enzyme) + 1
        py.subplot(2, 4, plotIndex)
        substrateDF = sampleDF[sampleDF["Enzyme"] == enzyme]
        py.plot(substrateDF["SubConcen"], substrateDF["Activity"])
        py.scatter(substrateDF["SubConcen"], substrateDF["Activity"])
        py.title("{0:}, {1:}, {2:}, {3:}".format(sample, vegetation,
                                                 precip, enzyme))
        py.xlabel("Substrate concentration (micromolar)")
        py.ylabel("Activity (micromole L^-1 g^-1 h^-1)")
    # Saving figures for data quality control purposes
    figName = figTitle + ".png"
    figPath = T0graphsFolder/figName
    py.savefig(figPath)
# %%
# Now, I'm going to begin working on T0 clear.
'''
T0ClearPath = enzymeFolderPath/enzymeFiles[0]
T0Clear = pd.read_csv(T0ClearPath, sep="\t", header=None, names=ogEnzymeCols)
'''
