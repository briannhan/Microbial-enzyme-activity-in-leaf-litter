# -*- coding: utf-8 -*-
"""
Created on Thu May 27 18:33:56 2021

@author: Brian Chung
This script creates boxplots of how Vmax and Km vary based on treatments and
time, and labels the boxplots with Tukey letters. Color schemes are meant to be
color-blind-friendly and are taken from the following link:
https://personal.sron.nl/~pault/#sec:qualitative

March 5, 2022 onwards
This script also creates boxplots for litter chemistry data
"""

import os
import pandas as pd
from matplotlib import pyplot as py
from matplotlib.patches import Patch
from pathlib import Path
from datetime import datetime
start = datetime.now()

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
'''I'm changing the treatments up. Originally timePoint was recorded as 0, 3,
5, and 6, but after briefly looking at a draft of a manuscript by Ashish, I'm
changing them so that they are, accordingly, 1, 2, 3, 4 to fall in line with
his manuscript. In addition, I'm shortening some of the treatment names
'''
renameTreatments = {"0": "1", "3": "2", "5": "3", "6": "4",  # timepoints only
                    "0 x Grassland": "1, Gr", "3 x Grassland": "2, Gr",
                    "5 x Grassland": "3, Gr", "6 x Grassland": "4, Gr",
                    "0 x CSS": "1, CSS", "3 x CSS": "2, CSS",
                    "5 x CSS": "3, CSS", "6 x CSS": "4, CSS",
                    # timePoint x Vegetation
                    "0 x Reduced": "1, D", "3 x Reduced": "2, D",
                    "5 x Reduced": "3, D", "6 x Reduced": "4, D",
                    "0 x Ambient": "1, A", "3 x Ambient": "2, A",
                    "5 x Ambient": "3, A", "6 x Ambient": "4, A",
                    # timePoint x Precipitation
                    "0 x Reduced x Grassland": "1, Gr, D",  # D = Drought
                    "3 x Reduced x Grassland": "2, Gr, D",
                    "5 x Reduced x Grassland": "3, Gr, D",
                    "6 x Reduced x Grassland": "4, Gr, D",
                    "0 x Ambient x Grassland": "1, Gr, A",  # A = Ambient
                    "3 x Ambient x Grassland": "2, Gr, A",
                    "5 x Ambient x Grassland": "3, Gr, A",
                    "6 x Ambient x Grassland": "4, Gr, A",
                    "0 x Reduced x CSS": "1, CSS, D",
                    "3 x Reduced x CSS": "2, CSS, D",
                    "5 x Reduced x CSS": "3, CSS, D",
                    "6 x Reduced x CSS": "4, CSS, D",
                    "0 x Ambient x CSS": "1, CSS, A",
                    "3 x Ambient x CSS": "2, CSS, A",
                    "5 x Ambient x CSS": "3, CSS, A",
                    "6 x Ambient x CSS": "4, CSS, A",
                    # three-way (PPO Vmax)
                    "CSS x Reduced": "CSS, D",
                    "CSS x Ambient": "CSS, A",
                    "Grassland x Reduced": "Gr, D",
                    "Grassland x Ambient": "Gr, A"}
# Vegetation x Precipitation
oldTreatments = renameTreatments.keys()
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
        newGroup = renameTreatments[oldGroup]
    annotationsAP_Km.loc[index, "groups"] = newGroup

# (3) Wrangling the parameters dataframe to create columns of Tukey groups
boolean1 = (parameters.Parameter == "Km") & (parameters.Enzyme == "AP")
AP_Km = parameters[boolean1]
paramColumns = AP_Km.columns.tolist()
splitAnnotateName = annotationsFileName.split(", ")
mainEorInteraction = splitAnnotateName[1]
indeVars = ["timePoint", "Vegetation", "Precipitation"]
# Making a Tukey Group column and then pivoting it to use its values as the
# new columns for the isolated parameter values
if "x" in mainEorInteraction:
    splittedFactors = mainEorInteraction.split(" x ")
    factor1 = splittedFactors[0]
    factor2 = splittedFactors[1]
    groupSeries = AP_Km[factor1] + " x " + AP_Km[factor2]
elif mainEorInteraction == "Three-way":
    groupSeries = (AP_Km["timePoint"] + " x " + AP_Km["Precipitation"] + " x "
                   + AP_Km["Vegetation"])
elif mainEorInteraction in indeVars:
    groupSeries = AP_Km[mainEorInteraction]
AP_Km["TukeyGroups"] = groupSeries
colsToRemove = []
for column in paramColumns:
    if column not in ["TukeyGroups", "value"]:
        colsToRemove.append(column)
AP_Km = AP_Km.drop(columns=colsToRemove)
AP_Km = AP_Km.pivot(columns="TukeyGroups", values="value")
AP_Km = AP_Km.rename(columns=renameTreatments)
TukeyGroups = annotationsAP_Km["groups"].tolist()
TukeyLabels = annotationsAP_Km["labels"].tolist()

# (4) Plotting boxplot for AP Km timePoint x Vegetation
dataToPlot = []
for group in TukeyGroups:
    dataColumn = AP_Km[group].dropna()
    dataToPlot.append(dataColumn)
AP_Kmcols = AP_Km.columns.tolist()
py.figure("AP, Km, timePoint x Vegetation", (20, 10))
bp = py.boxplot(dataToPlot, patch_artist=True, labels=TukeyGroups,
                showfliers=False)
y = r"Reaction products, $Log_{10}$ $K_m$ ($log_{10}$ $(\mu M)$)"
py.ylabel(y, fontfamily="serif", fontsize="xx-large", fontstyle="oblique")
py.yticks(fontsize="xx-large")
py.xlabel("Time, Vegetation combination", fontfamily="serif",
          fontsize="xx-large", fontstyle="oblique")
py.xticks(fontsize="xx-large")
py.title("AP, Km, Time, Vegetation combination", fontfamily="serif",
         fontsize="xx-large", fontstyle="oblique")
numOfTukeyGroups = len(bp["boxes"])
# Annotating Tukey labels
for i in range(numOfTukeyGroups):
    capIndex = (i*2) + 1
    topCap = bp["caps"][capIndex]
    xCap = topCap.get_xdata()
    yCap = topCap.get_ydata()
    xMiddle = xCap[0] + ((xCap[1] - xCap[0])/2)
    yVal = yCap[0]
    py.annotate(TukeyLabels[i], (xMiddle, yVal), (3, 4), "data",
                "offset pixels",
                fontsize="large", fontstyle="italic", fontweight="bold")
# Filling boxes in the boxplot & thickening the median line of each box
CSScolor = (34/255, 136/255, 51/255)  # green
grassColor = (204/255, 187/255, 68/255)  # yellow
white = (1, 1, 1)
# medianColor = (102/255, 204/255, 238/255)  # cyan
medianColor = (0, 0, 0)  # black
ambientHatch = "."
reducedHatch = "//"
for i in range(numOfTukeyGroups):
    series = dataToPlot[i]
    group = series.name
    patch = bp["boxes"][i]
    # Setting colors for vegetation type
    if "CSS" in group:
        faceColor = CSScolor
    elif "Gr" in group:
        faceColor = grassColor
    else:
        faceColor = white
    # Setting fill type for precipitation treatment
    if "A" in group:
        hatchType = ambientHatch
    elif "R" in group:
        hatchType = reducedHatch
    else:
        hatchType = None
    patch.set(facecolor=faceColor, hatch=hatchType)
    median = bp["medians"][i]
    median.set_linewidth(3)
    median.set_color(medianColor)
# Making a legend
CSSpatch = Patch(facecolor=CSScolor)
grassPatch = Patch(facecolor=grassColor)
ambientLegend = Patch(facecolor=white, hatch=ambientHatch)
reducedLegend = Patch(facecolor=white, hatch=reducedHatch)
if "Vegetation" in mainEorInteraction:
    vegPatches = [CSSpatch, grassPatch]
    vegLabels = ["Coastal sage scrub", "Grassland"]
else:
    vegPatches = []
    vegLabels = []
if "Precipitation" in mainEorInteraction:
    pptPatches = [ambientLegend, reducedLegend]
    pptLabels = ["Ambient", "Drought"]
else:
    pptPatches = []
    pptLabels = []
patches = vegPatches + pptPatches
legendLabels = vegLabels + pptLabels
py.legend(handles=patches, labels=legendLabels)
# py.savefig(APtukeyFolder/"AP, Km, timePoint x Vegetation",
#             bbox_inches="tight")


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
        if "annotated" in file and "worked out" not in file:
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
            newGroup = renameTreatments[oldGroup]
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
    # for index, row in paramsOI.iterrows():
    #     oldGroup = row["timePoint"]
    #     if oldGroup in oldTreatments:
    #         paramsOI.loc[index, "timePoint"] = renameTreatments[oldGroup]
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
    paramsOI = paramsOI.rename(columns=renameTreatments)
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
    paramsOIcols = paramsOI.columns.tolist()
    assert paramsOIcols == TukeyGroups
    for group in TukeyGroups:
        dataColumn = paramsOI[group].dropna()
        dataToPlot.append(dataColumn)
    fileNameSplit = fileName.split(", ")
    parameter = fileNameSplit[0]
    mainEorInteraction = fileNameSplit[1]
    figName = "{0}, {1}, {2}".format(enzyme, parameter, mainEorInteraction)
    py.figure(figName, (20, 10))
    bp = py.boxplot(dataToPlot, patch_artist=True, labels=TukeyGroups,
                    showfliers=False)
    if parameter == "Km":
        y = r"$Log_{10}$ $K_m$ ($log_{10}$ $(\mu M)$)"
        formattedParam = r"$K_{m}$"
    elif parameter == "Vmax":
        y = r"Enzyme amount, $Log_{10}$ $V_{max}$ ($log_{10}$ $(\mu mol/g/h)$)"
        formattedParam = r"$V_{max}$"
    py.ylabel(y, fontfamily="serif", fontsize="xx-large",
              fontstyle="oblique")
    py.yticks(fontsize="xx-large")

    if mainEorInteraction == "timePoint x Vegetation":
        xAxis = "Time, Vegetation combination"
    elif mainEorInteraction == "timePoint x Precipitation":
        xAxis = "Time, Precipitation combination"
    elif mainEorInteraction == "Vegetation x Precipitation":
        xAxis = "Vegetation, Precipitation combination"
    elif mainEorInteraction == "Three-way":
        xAxis = "Time & treatment combination"
    elif mainEorInteraction == "timePoint":
        xAxis = "Time"
    else:
        xAxis = mainEorInteraction
    py.xlabel(xAxis, fontfamily="serif", fontsize="xx-large",
              fontstyle="oblique")
    if mainEorInteraction != "Three-way":
        py.xticks(fontsize="xx-large")
    elif mainEorInteraction == "Three-way":
        py.xticks(fontsize="large")
    plotTitle = r"{0}, {1}, {2}".format(enzyme, formattedParam, xAxis)
    py.title(plotTitle, fontfamily="serif", fontsize="xx-large",
             fontstyle="oblique")

    # Annotating Tukey labels
    numOfTukeyGroups = len(TukeyGroups)
    for i in range(numOfTukeyGroups):
        capIndex = (i*2) + 1
        topCap = bp["caps"][capIndex]
        xCap = topCap.get_xdata()
        yCap = topCap.get_ydata()
        xMiddle = xCap[0] + ((xCap[1] - xCap[0])/2)
        yVal = yCap[0]
        py.annotate(TukeyLabels[i], (xMiddle, yVal), (5, 4), "data",
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
        elif "Gr" in group:
            faceColor = grassColor
        else:
            faceColor = white
        # Setting fill type for precipitation treatment
        if "A" in group:
            hatchType = ambientHatch
        elif "D" in group:
            hatchType = reducedHatch
        else:
            hatchType = None
        patch.set(facecolor=faceColor, hatch=hatchType)
        median = bp["medians"][i]
        median.set_linewidth(3)
        median.set_color(medianColor)

    # Making a legend
    all3 = "Three-way"
    if "Vegetation" in mainEorInteraction or all3 in mainEorInteraction:
        vegPatches = [CSSpatch, grassPatch]
        vegLabels = ["Coastal sage scrub", "Grassland"]
    else:
        vegPatches = []
        vegLabels = []
    if "Precipitation" in mainEorInteraction or all3 in mainEorInteraction:
        pptPatches = [ambientLegend, reducedLegend]
        pptLabels = ["Ambient", "Drought"]
    else:
        pptPatches = []
        pptLabels = []
    patches = vegPatches + pptPatches
    legendLabels = vegLabels + pptLabels
    if len(legendLabels) > 0 and len(patches) > 0:
        py.legend(handles=patches, labels=legendLabels)

    figPath = tukeyFolder/enzyme/figName
    py.savefig(figPath, bbox_inches="tight", pad_inches=0.04)
    return


def plotAllBoxPlots(enzyme):
    """
    Plot all the box plots associated with all the interactions and/or main
    effects associated with a particular enzyme.

    Parameters
    ----------
    enzyme : str
        Enzyme of interest. Will be used to choose the folder that contains the
        Tukey results for this enzyme as well as the parameter values specific
        to this enzyme.

    Returns
    -------
    None.

    """
    annotateFiles = annotationFiles(enzyme)
    numOfBoxPlots = len(annotateFiles)
    for i in range(numOfBoxPlots):
        plotBoxPlot(enzyme, annotateFiles[i])
    return


# %%
# Purpose: Making remaining boxplots for AP & for remaining enzymes

# (1) Making remaining boxplots for AP
# plotBoxPlot("AP", 'Km, Vegetation x Precipitation, groups, annotated.xlsx')
# plotBoxPlot("AP", 'Vmax, timePoint, groups, annotated.xlsx')

# (2) Making boxplots for all enzymes
# plotAllBoxPlots("BG")
# plotAllBoxPlots("BX")
# plotAllBoxPlots("CBH")
# plotAllBoxPlots("LAP")
# plotAllBoxPlots("NAG")
# plotAllBoxPlots("PPO")
# plotBoxPlot("PPO", "Vmax, Vegetation, groups, annotated.xlsx")
# %%
print(datetime.now() - start)
