# -*- coding: utf-8 -*-
"""
Created on Thu May 27 18:33:56 2021

@author: Brian Chung
This script creates boxplots of how Vmax and Km vary based on treatments and
time, and labels the boxplots with Tukey letters. Color schemes are meant to be
color-blind-friendly and are taken from the following link:
https://personal.sron.nl/~pault/#sec:qualitative
"""

import os
import pandas as pd
from matplotlib import pyplot as py
from matplotlib.patches import Patch
from pathlib import Path

# Finding paths that contain Tukey results
cwd = Path(os.getcwd())
cwdFilesDirs = os.listdir(cwd)
statsFolder = cwd/'Statistical analyses'
statsContents = os.listdir(statsFolder)
tukeyFolder = statsFolder/'Tukey posthoc'
tukeyContents = os.listdir(tukeyFolder)

# Loading in Michaelis-Menten parameters dataframe
activityFolder = cwd/'Enzyme activity data'
activityContents = os.listdir(activityFolder)
paramPath = activityFolder/'Parameters - log 10 transformed.xlsx'
parameters = pd.read_excel(paramPath).drop(columns=["Replicate",
                                                    "Transformation"])
parameters = parameters.rename(columns={"Precip": "Precipitation"})
parameters["timePoint"] = parameters["timePoint"].astype(str)

'''I'm changing the time points up. Originally they were recorded as 0, 3, 5,
and 6, but after briefly looking at a draft of a manuscript by Ashish, I'm
changing them so that they are, accordingly, 1, 2, 3, 4 to fall in line with
his manuscript.
'''
bool0_1 = parameters["timePoint"] == 0  # 0 is the o.g. timepoint, 1 is the new
parameters.loc[bool0_1, "timePoint"] = 1
bool3_2 = parameters["timePoint"] == 3  # 3 is the o.g. timepoint, 2 is the new
parameters.loc[bool3_2, "timePoint"] = 2
bool5_3 = parameters["timePoint"] == 5  # 5 is the o.g. timepoint, 3 is the new
parameters.loc[bool5_3, "timePoint"] = 3
bool6_4 = parameters["timePoint"] == 6  # 6 is the o.g. timepoint, 4 is the new
parameters.loc[bool6_4, "timePoint"] = 4
tpDict = {"0": "1", "3": "2", "5": "3", "6": "4",  # timepoints only
          "0 x Grassland": "1 x Grassland", "3 x Grassland": "2 x Grassland",
          "5 x Grassland": "3 x Grassland", "6 x Grassland": "4 x Grassland",
          "0 x CSS": "1 x CSS", "3 x CSS": "2 x CSS", "5 x CSS": "3 x CSS",
          "6 x CSS": "4 x CSS",  # timePoint x Vegetation
          "0 x Reduced": "1 x Reduced", "3 x Reduced": "2 x Reduced",
          "5 x Reduced": "3 x Reduced", "6 x Reduced": "4 x Reduced",
          "0 x Ambient": "1 x Ambient", "3 x Ambient": "2 x Ambient",
          "5 x Ambient": "3 x Ambient", "6 x Ambient": "4 x Ambient",
          # timePoint x Precipitation
          "0 x Reduced x Grassland": "1 x Grassland x Reduced",
          "3 x Reduced x Grassland": "2 x Grassland x Reduced",
          "5 x Reduced x Grassland": "3 x Grassland x Reduced",
          "6 x Reduced x Grassland": "4 x Grassland x Reduced",
          "0 x Ambient x Grassland": "1 x Grassland x Ambient",
          "3 x Ambient x Grassland": "2 x Grassland x Ambient",
          "5 x Ambient x Grassland": "3 x Grassland x Ambient",
          "6 x Ambient x Grassland": "4 x Grassland x Ambient",
          "0 x Reduced x CSS": "1 x CSS x Reduced",
          "3 x Reduced x CSS": "2 x CSS x Reduced",
          "5 x Reduced x CSS": "3 x CSS x Reduced",
          "6 x Reduced x CSS": "4 x CSS x Reduced",
          "0 x Ambient x CSS": "1 x CSS x Ambient",
          "3 x Ambient x CSS": "2 x CSS x Ambient",
          "5 x Ambient x CSS": "3 x CSS x Ambient",
          "6 x Ambient x CSS": "4 x CSS x Ambient"}  # three-way (PPO Vmax)
oldTreatments = tpDict.keys()
# This dict records the old time points & treatment combinations as the keys
# and the new time points & treatment combinations as the values
# %%
# Purpose: Plotting AP Vmax and/or Km possibly as a proof of concept to see
# if the process can be abstracted into a single function or a few functions

# (1) Isolating the files that contain the annotated groups
APtukeyFolder = tukeyFolder/"AP"
APtukeyContents = os.listdir(APtukeyFolder)
annotatedFiles = []
for file in APtukeyContents:
    if "annotated" in file:
        annotatedFiles.append(file)

# (2) Renaming treatment combinations/groups from the AP Km Tukey annotations
# Excel file
annotationsFileName = 'Km, timePoint x Vegetation, groups, annotated.xlsx'
annPath = APtukeyFolder/annotationsFileName
annotationsAP_Km = pd.read_excel(annPath).dropna(1)
annotationsAP_Km["groups"] = annotationsAP_Km["groups"].astype(str)
annotationsAP_Km = annotationsAP_Km.sort_values(by="groups")
for index, row in annotationsAP_Km.iterrows():
    oldGroup = row["groups"]
    if oldGroup not in oldTreatments:
        # If the annotation groups are either only vegetation, precipitation,
        # or vegetation x precipitation, then it exits this for loop
        break
    else:
        # If the annotation groups contain old time points, then they will be
        # renamed according to the dictionary above
        newGroup = tpDict[oldGroup]
    annotationsAP_Km.loc[index, "groups"] = newGroup

# (3) Wrangling the parameters dataframe to create columns of Tukey groups
boolean1 = (parameters.Parameter == "Km") & (parameters.Enzyme == "AP")
AP_Km = parameters[boolean1]
paramColumns = AP_Km.columns.tolist()
for index, row in AP_Km.iterrows():
    oldGroup = row["timePoint"]
    if oldGroup in oldTreatments:
        AP_Km.loc[index, "timePoint"] = tpDict[oldGroup]
splitAnnotateName = annotationsFileName.split(", ")
mainEorInteraction = splitAnnotateName[1]
if "x" in mainEorInteraction:
    splittedFactors = mainEorInteraction.split(" x ")
    factor1 = splittedFactors[0]
    factor2 = splittedFactors[1]
    groupSeries = AP_Km[factor1] + " x " + AP_Km[factor2]
elif mainEorInteraction == "Three-way":
    groupSeries = (AP_Km["timePoint"] + " x " + AP_Km["Precipitation"] + " x "
                   + AP_Km["Vegetation"])
else:
    groupSeries = AP_Km[mainEorInteraction]
AP_Km["TukeyGroups"] = groupSeries
colsToRemove = []
for column in paramColumns:
    if column not in ["TukeyGroups", "value"]:
        colsToRemove.append(column)
AP_Km = AP_Km.drop(columns=colsToRemove)
AP_Km = AP_Km.pivot(columns="TukeyGroups", values="value")
TukeyGroups = annotationsAP_Km["groups"].tolist()
TukeyLabels = annotationsAP_Km["labels"].tolist()

# (4) Plotting boxplot for AP Km timePoint x Vegetation
dataToPlot = []
for group in TukeyGroups:
    dataColumn = AP_Km[group].dropna()
    dataToPlot.append(dataColumn)
# py.figure("AP, Km, timePoint x Vegetation", (20, 14))
# bp = py.boxplot(dataToPlot, patch_artist=True, labels=TukeyGroups)
# py.ylabel("Log 10 Km (log10(micromolar))", fontfamily="serif",
#           fontsize="x-large", fontstyle="oblique")
# py.xlabel("Time, Vegetation combination", fontfamily="serif",
#           fontsize="x-large", fontstyle="oblique")
# numOfTukeyGroups = len(TukeyGroups)
# Annotating Tukey labels
# for i in range(numOfTukeyGroups):
#     capIndex = (i*2) + 1
#     topCap = bp["caps"][capIndex]
#     xCap = topCap.get_xdata()
#     yCap = topCap.get_ydata()
#     xMiddle = xCap[0] + ((xCap[1] - xCap[0])/2)
#     y = yCap[0]
#     py.annotate(TukeyLabels[i], (xMiddle, y), (3, 4), "data",
#                 "offset pixels",
#                 fontsize="large", fontstyle="italic", fontweight="bold")
# Filling boxes in the boxplot
CSScolor = (34/255, 136/255, 51/255)  # green
grassColor = (204/255, 187/255, 68/255)  # yellow
white = (1, 1, 1)
ambientHatch = "*"
reducedHatch = "//"
# for i in range(numOfTukeyGroups):
#     series = dataToPlot[i]
#     group = series.name
#     patch = bp["boxes"][i]
#     # Setting colors for vegetation type
#     if "CSS" in group:
#         faceColor = CSScolor
#     elif "Grassland" in group:
#         faceColor = grassColor
#     else:
#         faceColor = white
#     # Setting fill type for precipitation treatment
#     if "Ambient" in group:
#         hatchType = ambientHatch
#     elif "Reduced" in group:
#         hatchType = reducedHatch
#     else:
#         hatchType = None
#     patch.set(facecolor=faceColor, hatch=hatchType)
# Making a legend
CSSpatch = Patch(facecolor=CSScolor)
grassPatch = Patch(facecolor=grassColor)
ambientLegend = Patch(facecolor=white, hatch=ambientHatch)
reducedLegend = Patch(facecolor=white, hatch=reducedHatch)
# if "Vegetation" in mainEorInteraction:
#     vegPatches = [CSSpatch, grassPatch]
#     vegLabels = ["Coastal sage scrub", "Grassland"]
# else:
#     vegPatches = []
#     vegLabels = []
# if "Precipitation" in mainEorInteraction:
#     pptPatches = [ambientLegend, reducedLegend]
#     pptLabels = ["Ambient", "Reduced"]
# else:
#     pptPatches = []
#     pptLabels = []
# patches = vegPatches + pptPatches
# legendLabels = vegLabels + pptLabels
# py.legend(handles=patches, labels=legendLabels)
# py.savefig(APtukeyFolder/"AP, Km, timePoint x Vegetation")


# (5) Abstracting the process into a few functions
def annotationFiles(enzyme):
    """
    Looks into a folder that contains Tukey results for a particular enzyme
    and returns a list of Excel files with Tukey labels from that folder.

    Parameters
    ----------
    enzyme : str
        Enzyme of interest. This argument is used to identify the enzyme whose
        folder will be looked into.

    Returns
    -------
    List of strings where each string is an Excel file of Tukey labels from a
    particular main effect or interaction.

    """
    ezTukeyFolder = tukeyFolder/enzyme
    ezTukeyContents = os.listdir(ezTukeyFolder)
    annotatedFiles = []
    for file in ezTukeyContents:
        if "annotated" in file:
            annotatedFiles.append(file)
    return annotatedFiles


def annotationFileRename(enzyme, fileName):
    """
    Renames the treatment combinations that contains timePoint in a Tukey
    annotations file and returns a dataframe of the file.

    Parameters
    ----------
    enzyme : str
        Enzyme of interest. Used to look into the folder that contains Tukey
        results of a particular enzyme.
    fileName : str
        Name of the Tukey annotations file whose treatment combinations contain
        timePoint and are intended to be renamed.

    Returns
    -------
    Pandas dataframe of the Tukey annotations file with the treatment
    combinations renamed and sorted in ascending order.

    """
    ezTukeyFolder = tukeyFolder/enzyme
    annPath = ezTukeyFolder/fileName
    annotations = pd.read_excel(annPath).dropna(1)
    annotations["groups"] = annotations["groups"].astype(str)
    annotations = annotations.sort_values(by="groups")
    for index, row in annotations.iterrows():
        oldGroup = row["groups"]
        if oldGroup not in oldTreatments:
            # If the annotation groups are either only vegetation,
            # precipitation, or vegetation x precipitation, then it exits this
            # for loop
            break
        else:
            # If the annotation groups contain old time points, then they will
            # be renamed according to the dictionary above
            newGroup = tpDict[oldGroup]
        annotations.loc[index, "groups"] = newGroup
    return annotations


def TukeyGroupCols(enzyme, fileName):
    """
    Creates Tukey group columns out of the original independent variables in
    the parameter dataframe.

    Parameters
    ----------
    enzyme : str
        Enzyme of interest. Will be used to select the folder that contains the
        Tukey results for that enzyme as well as to select the parameter values
        for that enzyme.
    fileName : str
        Name of the annotation file. Will be used to select the parameter of
        interest as well as the main effect or interactions of interest to make
        Tukey group columns.

    Returns
    -------
    None.

    """
    splitAnnotateName = fileName.split(", ")
    parameter = splitAnnotateName[0]
    mainEorInteraction = splitAnnotateName[1]
    condi = (parameters.Parameter == parameter) & (parameters.Enzyme == enzyme)
    paramsOI = parameters[condi]
    paramOIogCols = paramsOI.columns.tolist()
    for index, row in paramsOI.iterrows():
        oldGroup = row["timePoint"]
        if oldGroup in oldTreatments:
            paramsOI.loc[index, "timePoint"] = tpDict[oldGroup]
    if "x" in mainEorInteraction:
        splittedFactors = mainEorInteraction.split(" x ")
        factor1 = splittedFactors[0]
        factor2 = splittedFactors[1]
        groupSeries = paramsOI[factor1] + " x " + paramsOI[factor2]
    elif mainEorInteraction == "Three-way":
        groupSeries = (paramsOI["timePoint"] + " x "
                       + paramsOI["Precipitation"] + " x "
                       + paramsOI["Vegetation"])
    else:
        groupSeries = paramsOI[mainEorInteraction]
    paramsOI["TukeyGroups"] = groupSeries
    colsToRemove = []
    for column in paramOIogCols:
        if column not in ["TukeyGroups", "value"]:
            colsToRemove.append(column)
    paramsOI = paramsOI.drop(columns=colsToRemove)
    paramsOI = paramsOI.pivot(columns="TukeyGroups", values="value")
    return paramsOI


def plotBoxPlot(enzyme, fileName):
    """
    Makes a boxplot of a parameter grouped by either an independent variable or
    interactions of at least 2 independent variables. This function calls on
    annotationFileRename() and TukeyGroupCols(), so these 2 functions do not
    need to be called outside of this function.

    Parameters
    ----------
    enzyme : str
        Enzyme of interest. Will be used to select the folder that contains
        Tukey results for this particular enzyme, to select the parameter
        values specific to this enzyme, and to name the boxplot.
    fileName : str
        Name of the annotation file. If the groups contain timePoint as an
        independent variable, then the groups will be renamed. In addition,
        this argument will be used to choose the parameter of interest and the
        independent variable(s) to make Tukey group columns.

    Returns
    -------
    None.

    """
    annotations = annotationFileRename(enzyme, fileName)
    paramsOI = TukeyGroupCols(enzyme, fileName)
    TukeyGroups = annotations["groups"].tolist()
    TukeyLabels = annotations["labels"].tolist()
    dataToPlot = []
    for group in TukeyGroups:
        dataColumn = paramsOI[group].dropna()
        dataToPlot.append(dataColumn)
    fileNameSplit = fileName.split(", ")
    parameter = fileNameSplit[0]
    mainEorInteraction = fileNameSplit[1]
    figName = "{0}, {1}, {2}".format(enzyme, parameter, mainEorInteraction)
    py.figure(figName, (20, 14))
    bp = py.boxplot(dataToPlot, patch_artist=True, labels=TukeyGroups)
    if parameter == "Km":
        yAxis = "Log 10 Km (log10(micromolar))"
    elif parameter == "Vmax":
        yAxis = "Log 10 Vmax (log10(micromole/g/s))"
    py.ylabel(yAxis, fontfamily="serif", fontsize="x-large",
              fontstyle="oblique")

    if mainEorInteraction == "timePoint x Vegetation":
        xAxis = "Time, Vegetation combination"
    elif mainEorInteraction == "timePoint x Precipitation":
        xAxis = "Time, Precipitation combination"
    elif mainEorInteraction == "Vegetation x Precipitation":
        xAxis = "Vegetation, Precipitation combination"
    elif mainEorInteraction == "Three-way":
        xAxis = "Three-way combination"
    py.xlabel(xAxis, fontfamily="serif", fontsize="x-large",
              fontstyle="oblique")

    # Annotating Tukey labels
    numOfTukeyGroups = len(TukeyGroups)
    for i in range(numOfTukeyGroups):
        capIndex = (i*2) + 1
        topCap = bp["caps"][capIndex]
        xCap = topCap.get_xdata()
        yCap = topCap.get_ydata()
        xMiddle = xCap[0] + ((xCap[1] - xCap[0])/2)
        y = yCap[0]
        py.annotate(TukeyLabels[i], (xMiddle, y), (3, 4), "data",
                    "offset pixels", fontsize="large", fontstyle="oblique",
                    fontweight="bold")

    # Filling boxes in the boxplot
    for i in range(numOfTukeyGroups):
        series = dataToPlot[i]
        group = series.name
        patch = bp["boxes"][i]
        # Setting colors for vegetation type
        if "CSS" in group:
            faceColor = CSScolor
        elif "Grassland" in group:
            faceColor = grassColor
        else:
            faceColor = white
        # Setting fill type for precipitation treatment
        if "Ambient" in group:
            hatchType = ambientHatch
        elif "Reduced" in group:
            hatchType = reducedHatch
        else:
            hatchType = None
        patch.set(facecolor=faceColor, hatch=hatchType)

    # Making a legend
    if "Vegetation" in mainEorInteraction:
        vegPatches = [CSSpatch, grassPatch]
        vegLabels = ["Coastal sage scrub", "Grassland"]
    else:
        vegPatches = []
        vegLabels = []
    if "Precipitation" in mainEorInteraction:
        pptPatches = [ambientLegend, reducedLegend]
        pptLabels = ["Ambient", "Reduced"]
    else:
        pptPatches = []
        pptLabels = []
    patches = vegPatches + pptPatches
    legendLabels = vegLabels + pptLabels
    if len(legendLabels) > 0 and len(patches) > 0:
        py.legend(handles=patches, labels=legendLabels)

    figPath = tukeyFolder/enzyme/figName
    py.savefig(figPath)
    return
